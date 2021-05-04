using System;
using System.Collections.Generic;
using System.Linq;
using NMath;
using System.IO;
using Microsoft.CSharp;
using System.CodeDom.Compiler;
using ConsoleExtensions;

namespace FEM
{
    class FirstBoundaryCondition
    {
        public Func<double, double, double> function { get; private set; }
        public List<int> points { get; private set; }
        public FirstBoundaryCondition(Func<double, double, double> function, List<int> points)
        {
            this.function = function;
            this.points = points;
        }
        public bool containsPoint(int point)
        {
            return points.Contains(point);
        }
    }
    class SecondBoundaryCondition
    {
        public Func<double, double, double> function { get; private set; }
        public List<int> points { get; private set; }
        public SecondBoundaryCondition(Func<double, double, double> function, List<int> points)
        {
            this.function = function;
            this.points = points;
        }
        public bool containsEdge(int pointA, int pointB)
        {
            int size = points.Count - 1;
            for (int i = 0; i < size; i++)
                if (points[i] == pointA && points[i + 1] == pointB || points[i] == pointB && points[i + 1] == pointB)
                    return true;
            return false;
        }
    }
    class Mesh
    {
        public Vector2[] vertexes { get; private set; }
        public int[][] polygons { get; private set; }
        public static Mesh Parse(string filePath)
        {
            if (!File.Exists(filePath))
                throw new FileNotFoundException();
            Mesh mesh = new Mesh();
            List<Vector2> vertexes = new List<Vector2>();
            List<int[]> polygons = new List<int[]>();
            StreamReader reader = new StreamReader(File.OpenRead(filePath));
            while (!reader.EndOfStream)
            {
                string[] data = reader.ReadLine().Split(' ');
                if (data.Length == 0)
                    continue;
                switch (data[0])
                {
                    case "v":
                        vertexes.Add(new Vector2(double.Parse(data[1].Replace('.', ',')), double.Parse(data[2].Replace('.', ','))));
                        break;
                    case "f":
                        polygons.Add(new int[] { int.Parse(data[1]), int.Parse(data[2]), int.Parse(data[3]) });
                        break;
                }
            }
            mesh.vertexes = vertexes.ToArray();
            mesh.polygons = polygons.ToArray();
            return mesh;
        }
    }
    class Solution
    {
        public Mesh mesh;
        public double[] coeffs;
        public List<KeyValuePair<int, int>>[] adjacencyList;
        public double get(double X, double Y)
        {
            bool isInside(int polygonIndex)
            {
                Vector2[] poly = mesh.polygons[polygonIndex].Select(index => mesh.vertexes[index]).ToArray();
                double Bx, By, Cx, Cy;
                Bx = poly[1].X - poly[0].X;
                By = poly[1].Y - poly[0].Y;
                Cx = poly[2].X - poly[0].X;
                Cy = poly[2].Y - poly[0].Y;
                double m = ((X - poly[0].X) * By - (Y - poly[0].Y) * Bx) / (Cx * By - Cy * Bx);
                if (m >= 0 && m <= 1)
                {
                    double l = ((X - poly[0].X) * Cy - (Y - poly[0].Y) * Cx) / (Bx * Cy - By * Cx);
                    if (l >= 0 && m + l <= 1)
                        return true;
                }
                return false;
            }
            int getAdditionalVertexBetween(int i, int j)
            {
                if (j > i)
                {
                    i = i + j;
                    j = i - j;
                    i = i - j;
                }
                foreach (KeyValuePair<int, int> pair in adjacencyList[i])
                    if (pair.Key == j)
                        return pair.Value;
                return -1;
            }
            for (int i = 0; i < mesh.polygons.Length; i++)
            {
                if (isInside(i))
                {
                    Vector2 A = mesh.vertexes[mesh.polygons[i][0]];
                    Vector2 B = mesh.vertexes[mesh.polygons[i][1]];
                    Vector2 C = mesh.vertexes[mesh.polygons[i][2]];
                    double detD = (B.X * C.Y - C.X * B.Y) - (A.X * C.Y - C.X * A.Y) + (A.X * B.Y - B.X * A.Y);
                    double[,] alphaMatrix = new double[3, 3];
                    alphaMatrix[0, 0] = (B.X * C.Y - C.X * B.Y) / detD;
                    alphaMatrix[1, 0] = -(A.X * C.Y - C.X * A.Y) / detD;
                    alphaMatrix[2, 0] = (A.X * B.Y - B.X * A.Y) / detD;
                    alphaMatrix[0, 1] = -(C.Y - B.Y) / detD;
                    alphaMatrix[1, 1] = (C.Y - A.Y) / detD;
                    alphaMatrix[2, 1] = -(B.Y - A.Y) / detD;
                    alphaMatrix[0, 2] = (C.X - B.X) / detD;
                    alphaMatrix[1, 2] = -(C.X - A.X) / detD;
                    alphaMatrix[2, 2] = (B.X - A.X) / detD;

                    double L1, L2, L3;
                    L1 = alphaMatrix[0, 0] + alphaMatrix[0, 1] * X + alphaMatrix[0, 2] * Y;
                    L2 = alphaMatrix[1, 0] + alphaMatrix[1, 1] * X + alphaMatrix[1, 2] * Y;
                    L3 = alphaMatrix[2, 0] + alphaMatrix[2, 1] * X + alphaMatrix[2, 2] * Y;

                    double psi1, psi2, psi3, psi4, psi5, psi6;
                    psi1 = 2.0 * L1 * L1 - L1;
                    psi2 = 2.0 * L2 * L2 - L2;
                    psi3 = 2.0 * L3 * L3 - L3;
                    psi4 = 4.0 * L1 * L2;
                    psi5 = 4.0 * L2 * L3;
                    psi6 = 4.0 * L1 * L3;

                    return coeffs[mesh.polygons[i][0]] * psi1 +
                           coeffs[mesh.polygons[i][1]] * psi2 +
                           coeffs[mesh.polygons[i][2]] * psi3 +
                           coeffs[getAdditionalVertexBetween(mesh.polygons[i][0], mesh.polygons[i][1])] * psi4 +
                           coeffs[getAdditionalVertexBetween(mesh.polygons[i][1], mesh.polygons[i][2])] * psi5 +
                           coeffs[getAdditionalVertexBetween(mesh.polygons[i][2], mesh.polygons[i][0])] * psi6;
                }
            }
            throw new ArgumentOutOfRangeException();
        }
    }
    class Problem
    {
        public Mesh mesh;
        public Func<double, double, double> f;
        public List<FirstBoundaryCondition> firstBoundaryConditions = new List<FirstBoundaryCondition>();
        public List<SecondBoundaryCondition> secondBoundaryConditions = new List<SecondBoundaryCondition>();
        public double[] lambda;
        public double[] gamma;

        private List<KeyValuePair<int, int>>[] adjacencyList;
        public Solution solve()
        {
            int additionalVertexesCount = mesh.vertexes.Length + mesh.polygons.Length - 1;

            int size = mesh.vertexes.Length + additionalVertexesCount;
            int[] ig = new int[size + 1];
            int[] jg = new int[mesh.polygons.Length * 6 + additionalVertexesCount * 3];
            double[] di = new double[size];
            double[] gg = new double[mesh.polygons.Length * 6 + additionalVertexesCount * 3];
            MatrixSymm_Sparse globalMatrix = new MatrixSymm_Sparse(size, ig, jg, di, gg);
            Vector globalVector = new Vector(size);
            double[] globalVec = globalVector.values;

            adjacencyList = new List<KeyValuePair<int, int>>[mesh.vertexes.Length];
            for (int i = 0; i < mesh.vertexes.Length; i++)
                adjacencyList[i] = new List<KeyValuePair<int, int>>();

            int getAdditionalVertexBetween(int i, int j)
            {
                if (j > i)
                {
                    i = i + j;
                    j = i - j;
                    i = i - j;
                }
                foreach (KeyValuePair<int, int> pair in adjacencyList[i])
                    if (pair.Key == j)
                        return pair.Value;
                return -1;
            }

            int additionalVertexIndex = mesh.vertexes.Length;
            for (int i = 0; i < mesh.vertexes.Length; i++)
            {
                ig[i + 1] = ig[i];
                for (int j = 0; j < i; j++)
                {
                    for (int k = 0; k < mesh.polygons.Length; k++)
                    {
                        int[] curPoly = mesh.polygons[k];
                        for (int w = 0; w < 3; w++)
                            if (curPoly[w] == i && curPoly[w == 2 ? 0 : w + 1] == j || curPoly[w] == j && curPoly[w == 2 ? 0 : w + 1] == i)
                            {
                                adjacencyList[i].Add(new KeyValuePair<int, int>(j, additionalVertexIndex));
                                additionalVertexIndex++;
                                jg[ig[i + 1]] = j;
                                ig[i + 1]++;

                                k = mesh.polygons.Length;
                                break;
                            }
                    }
                }
            }
            for (int i = mesh.vertexes.Length; i < size; i++)
            {
                ig[i + 1] = ig[i];
                for (int j = 0; j < i; j++)
                {
                    foreach (int[] polygon in mesh.polygons)
                    {
                        if (j < mesh.vertexes.Length)
                        {
                            if (polygon.Contains(j) && (getAdditionalVertexBetween(polygon[0], polygon[1]) == i ||
                                                        getAdditionalVertexBetween(polygon[1], polygon[2]) == i ||
                                                        getAdditionalVertexBetween(polygon[2], polygon[0]) == i))
                            {
                                jg[ig[i + 1]] = j;
                                ig[i + 1]++;
                                break;
                            }
                        }
                        else
                        {
                            int a, b, c;
                            a = getAdditionalVertexBetween(polygon[0], polygon[1]);
                            b = getAdditionalVertexBetween(polygon[1], polygon[2]);
                            c = getAdditionalVertexBetween(polygon[2], polygon[0]);
                            if ((a == i || b == i || c == i) && (a == j || b == j || c == j))
                            {
                                jg[ig[i + 1]] = j;
                                ig[i + 1]++;
                                break;
                            }
                        }
                    }
                }
            }

            double[][] localMatrix = new double[6][];
            for (int i = 0; i < 6; i++)
                localMatrix[i] = new double[i + 1];
            double[] localVector = new double[6];
            double[] fValues = new double[6];
            double[] conditionValues = new double[6];
            double[,] alphaMatrix = new double[3, 3];
            double detD, absdetD;
            Vector2 A, B, C;
            for (int p = 0; p < mesh.polygons.Length; p++)
            {
                A = mesh.vertexes[mesh.polygons[p][0]];
                B = mesh.vertexes[mesh.polygons[p][1]];
                C = mesh.vertexes[mesh.polygons[p][2]];
                detD = (B.X * C.Y - C.X * B.Y) - (A.X * C.Y - C.X * A.Y) + (A.X * B.Y - B.X * A.Y);
                absdetD = Math.Abs(detD);
                alphaMatrix[0, 0] = (B.X * C.Y - C.X * B.Y) / detD;
                alphaMatrix[1, 0] = -(A.X * C.Y - C.X * A.Y) / detD;
                alphaMatrix[2, 0] = (A.X * B.Y - B.X * A.Y) / detD;
                alphaMatrix[0, 1] = -(C.Y - B.Y) / detD;
                alphaMatrix[1, 1] = (C.Y - A.Y) / detD;
                alphaMatrix[2, 1] = -(B.Y - A.Y) / detD;
                alphaMatrix[0, 2] = (C.X - B.X) / detD;
                alphaMatrix[1, 2] = -(C.X - A.X) / detD;
                alphaMatrix[2, 2] = (B.X - A.X) / detD;

                // mass matrix
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j <= i; j++)
                    {
                        localMatrix[i][j] = (i == j ? 1.0 / 60.0 : -1.0 / 360.0) * absdetD;
                        localMatrix[i + 3][j + 3] = (i == j ? 4.0 / 45.0 : 2.0 / 45.0) * absdetD;
                    }
                    localMatrix[3][i] = (i != 2 ? 0 : -1.0 / 90.0) * absdetD;
                    localMatrix[4][i] = (i != 0 ? 0 : -1.0 / 90.0) * absdetD;
                    localMatrix[5][i] = (i != 1 ? 0 : -1.0 / 90.0) * absdetD;
                }

                fValues[0] = f(A.X, A.Y);
                fValues[1] = f(B.X, B.Y);
                fValues[2] = f(C.X, C.Y);
                fValues[3] = f((A.X + B.X) / 2.0, (A.Y + B.Y) / 2.0);
                fValues[4] = f((B.X + C.X) / 2.0, (B.Y + C.Y) / 2.0);
                fValues[5] = f((A.X + C.X) / 2.0, (A.Y + C.Y) / 2.0);
                // vector
                for (int i = 0; i < 6; i++)
                    localVector[i] = 0;
                for (int i = 0; i < 6; i++)
                {
                    localVector[i] += localMatrix[i][i] * fValues[i];
                    for (int j = 0; j < i; j++)
                    {
                        localVector[i] += localMatrix[i][j] * fValues[j];
                        localVector[j] += localMatrix[i][j] * fValues[i];

                        localMatrix[i][j] *= gamma[p];
                    }
                    localMatrix[i][i] *= gamma[p];
                }

                foreach (SecondBoundaryCondition condition in secondBoundaryConditions)
                {
                    conditionValues[0] = condition.function(A.X, A.Y);
                    conditionValues[1] = condition.function(B.X, B.Y);
                    conditionValues[2] = condition.function(C.X, C.Y);
                    conditionValues[3] = condition.function((A.X + B.X) / 2.0, (A.Y + B.Y) / 2.0);
                    conditionValues[4] = condition.function((B.X + C.X) / 2.0, (B.Y + C.Y) / 2.0);
                    conditionValues[5] = condition.function((A.X + C.X) / 2.0, (A.Y + C.Y) / 2.0);
                    int[] polygon = mesh.polygons[p];
                    int nextPointIndex;
                    for (int i = 0; i < 3; i++)
                        if (condition.containsEdge(polygon[0], polygon[1]))
                        {
                            nextPointIndex = i == 2 ? 0 : i + 1;
                            localVector[i] += (4.0 * conditionValues[i] +
                                               2.0 * conditionValues[i + 3] -
                                               conditionValues[nextPointIndex]) / 30.0;
                            localVector[nextPointIndex] += (4.0 * conditionValues[nextPointIndex] +
                                               2.0 * conditionValues[i + 3] -
                                               conditionValues[i]) / 30.0;
                            localVector[i + 3] += (2.0 * conditionValues[i] +
                                               16.0 * conditionValues[i + 3] +
                                               2.0 * conditionValues[nextPointIndex]) / 30.0;
                        }
                }

                // stiffness matrix
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j <= i; j++)
                        localMatrix[i][j] += lambda[p] * (i == j ? 1.0 / 2.0 : -1.0 / 6.0) * (alphaMatrix[i, 1] * alphaMatrix[j, 1] + alphaMatrix[i, 2] * alphaMatrix[j, 2]) * absdetD;
                    localMatrix[3][i] += lambda[p] * (((i == 1 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[0, 1] + alphaMatrix[i, 2] * alphaMatrix[0, 2])) +
                                                      ((i == 0 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[1, 1] + alphaMatrix[i, 2] * alphaMatrix[1, 2]))) * absdetD;
                    localMatrix[4][i] += lambda[p] * (((i == 2 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[1, 1] + alphaMatrix[i, 2] * alphaMatrix[1, 2])) +
                                                      ((i == 1 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[2, 1] + alphaMatrix[i, 2] * alphaMatrix[2, 2]))) * absdetD;
                    localMatrix[5][i] += lambda[p] * (((i == 2 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[0, 1] + alphaMatrix[i, 2] * alphaMatrix[0, 2])) +
                                                      ((i == 0 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[2, 1] + alphaMatrix[i, 2] * alphaMatrix[2, 2]))) * absdetD;
                }
                localMatrix[3][3] += 4.0 / 3.0 * lambda[p] * (alphaMatrix[0, 1] * alphaMatrix[0, 1] + alphaMatrix[0, 2] * alphaMatrix[0, 2] +
                                                              alphaMatrix[0, 1] * alphaMatrix[1, 1] + alphaMatrix[0, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[1, 1] + alphaMatrix[1, 2] * alphaMatrix[1, 2]) * absdetD;
                localMatrix[4][4] += 4.0 / 3.0 * lambda[p] * (alphaMatrix[1, 1] * alphaMatrix[1, 1] + alphaMatrix[1, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[2, 1] + alphaMatrix[1, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[2, 1] * alphaMatrix[2, 1] + alphaMatrix[2, 2] * alphaMatrix[2, 2]) * absdetD;
                localMatrix[5][5] += 4.0 / 3.0 * lambda[p] * (alphaMatrix[2, 1] * alphaMatrix[2, 1] + alphaMatrix[2, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[2, 1] * alphaMatrix[0, 1] + alphaMatrix[2, 2] * alphaMatrix[0, 2] +
                                                              alphaMatrix[0, 1] * alphaMatrix[0, 1] + alphaMatrix[0, 2] * alphaMatrix[0, 2]) * absdetD;
                localMatrix[4][3] += 2.0 / 3.0 * lambda[p] * (alphaMatrix[0, 1] * alphaMatrix[1, 1] + alphaMatrix[0, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[2, 1] + alphaMatrix[1, 2] * alphaMatrix[2, 2] +
                                                              2.0 * alphaMatrix[0, 1] * alphaMatrix[2, 1] + 2.0 * alphaMatrix[0, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[1, 1] + alphaMatrix[1, 2] * alphaMatrix[1, 2]) * absdetD;
                localMatrix[5][3] += 2.0 / 3.0 * lambda[p] * (alphaMatrix[0, 1] * alphaMatrix[1, 1] + alphaMatrix[0, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[0, 1] * alphaMatrix[2, 1] + alphaMatrix[0, 2] * alphaMatrix[2, 2] +
                                                              2.0 * alphaMatrix[1, 1] * alphaMatrix[2, 1] + 2.0 * alphaMatrix[1, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[0, 1] * alphaMatrix[0, 1] + alphaMatrix[0, 2] * alphaMatrix[0, 2]) * absdetD;
                localMatrix[5][4] += 2.0 / 3.0 * lambda[p] * (alphaMatrix[0, 1] * alphaMatrix[2, 1] + alphaMatrix[0, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[2, 1] + alphaMatrix[1, 2] * alphaMatrix[2, 2] +
                                                              2.0 * alphaMatrix[0, 1] * alphaMatrix[1, 1] + 2.0 * alphaMatrix[0, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[2, 1] * alphaMatrix[2, 1] + alphaMatrix[2, 2] * alphaMatrix[2, 2]) * absdetD;

                void addToGlobalMatrix(int row, int column, double value)
                {
                    if (row == column)
                    {
                        di[row] += value;
                        return;
                    }
                    if (row < column)
                    {
                        row = row + column;
                        column = row - column;
                        row = row - column;
                    }
                    for (int k = ig[row]; k < ig[row + 1]; k++)
                        if (jg[k] == column)
                        {
                            gg[k] += value;
                            return;
                        }
                }

                int[] globalIndexes = new int[6];
                for (int i = 0; i < 3; i++)
                {
                    globalIndexes[i] = mesh.polygons[p][i];
                    globalIndexes[i + 3] = getAdditionalVertexBetween(mesh.polygons[p][i], mesh.polygons[p][i == 2 ? 0 : i + 1]);
                }    
                for (int i = 0; i < 3; i++)
                {
                    globalVec[globalIndexes[i]] += localVector[i];
                    globalVec[globalIndexes[i + 3]] += localVector[i + 3];
                    for (int j = 0; j < 3; j++)
                    {
                        if (j <= i)
                        {
                            addToGlobalMatrix(globalIndexes[i], globalIndexes[j], localMatrix[i][j]);
                            addToGlobalMatrix(globalIndexes[i + 3], globalIndexes[j + 3], localMatrix[i + 3][j + 3]);
                        }
                        addToGlobalMatrix(globalIndexes[i + 3], globalIndexes[j], localMatrix[i + 3][j]);
                    }
                }
            }

            void setFirstCondition(int index, double value)
            {
                globalVec[index] = value;
                for (int i = ig[index]; i < ig[index + 1] && jg[i] < index; i++)
                {
                    globalVec[jg[i]] -= globalVec[index] * gg[i];
                    gg[i] = 0;
                }
                for (int i = index + 1; i < size; i++)
                    for (int j = ig[i]; j < ig[i + 1] && jg[j] <= index; j++)
                        if (jg[j] == index)
                        {
                            globalVec[i] -= globalVec[index] * gg[j];
                            gg[j] = 0;
                            break;
                        }
                di[index] = 1.0;
            }

            foreach (FirstBoundaryCondition condition in firstBoundaryConditions)
            {
                int count = condition.points.Count;
                int i;
                int pointIndex;
                int nextPointIndex;
                nextPointIndex = condition.points[0];
                for (i = 0; i < count - 1; i++)
                {
                    pointIndex = nextPointIndex;
                    nextPointIndex = condition.points[i + 1];
                    setFirstCondition(pointIndex, condition.function(mesh.vertexes[pointIndex].X, mesh.vertexes[pointIndex].Y));
                    setFirstCondition(getAdditionalVertexBetween(pointIndex, nextPointIndex), 
                        condition.function((mesh.vertexes[pointIndex].X + mesh.vertexes[nextPointIndex].X) / 2.0, 
                                           (mesh.vertexes[pointIndex].Y + mesh.vertexes[nextPointIndex].Y) / 2.0));
                }
                setFirstCondition(condition.points[i], condition.function(mesh.vertexes[condition.points[i]].X, mesh.vertexes[condition.points[i]].Y));
            }

            Vector coeffs = new Vector(size);
            SLAE.Solve_CGM_LLt(globalMatrix, globalVector, coeffs, out _);
            Solution solution = new Solution();
            solution.coeffs = coeffs.values;
            solution.adjacencyList = adjacencyList;
            solution.mesh = mesh;
            return solution;
        }
    }
    class Program
    {
        static void Main(string[] args)
        {
            // -div(lambda*grad(u)) + gamma*u = f
            Problem p = new Problem();
            p.mesh = Mesh.Parse("mesh.txt");
            string[] materials = File.ReadAllLines("materials.txt");
            p.lambda = materials[0].Split(' ').Select(value => double.Parse(value)).ToArray();
            p.gamma = materials[1].Split(' ').Select(value => double.Parse(value)).ToArray();

            p.f = ParseExpression(File.ReadAllText("function.txt"));
            string[] lines = File.ReadAllLines("first boundary.txt");
            for (int i = 0; i < lines.Length; i += 2)
                p.firstBoundaryConditions.Add(new FirstBoundaryCondition(ParseExpression(lines[i]), lines[i + 1].Split(' ').Select(value => int.Parse(value)).ToList()));
            lines = File.ReadAllLines("second boundary.txt");
            for (int i = 0; i < lines.Length; i += 2)
                p.secondBoundaryConditions.Add(new SecondBoundaryCondition(ParseExpression(lines[i]), lines[i + 1].Split(' ').Select(value => int.Parse(value)).ToList()));
            
            Solution solution = p.solve();

            int choice;
            while ((choice = ConsoleMenu.Show("Select option:", "get value at point", "get grid of values", "get vector of coeffs", "exit")) != 3)
            {
                switch (choice)
                {
                    case 0:
                        {
                            Console.WriteLine("Write point X and Y:");
                            string[] input = Console.ReadLine().Replace('.', ',').Split(' ');
                            try
                            {
                                Console.WriteLine(solution.get(double.Parse(input[0]), double.Parse(input[1])).ToString("E5"));
                            }
                            catch
                            {
                                Console.WriteLine("Point does not belong to the mesh area");
                            }

                            Console.WriteLine("Press any key to get back");
                            Console.ReadKey();
                            break;
                        }
                    case 1:
                        {
                            Console.WriteLine("Write start X, start Y, delta X, delta Y, amount of points on X and amount of points on Y:");
                            string[] input = Console.ReadLine().Replace('.', ',').Split(' ');
                            double startX, startY, deltaX, deltaY;
                            int countX, countY;
                            startX = double.Parse(input[0]);
                            startY = double.Parse(input[1]);
                            deltaX = double.Parse(input[2]);
                            deltaY = double.Parse(input[3]);
                            countX = int.Parse(input[4]);
                            countY = int.Parse(input[5]);
                            for (int i = 0; i < countY; i++)
                            {
                                for (int j = 0; j < countX; j++)
                                    try
                                    {
                                        Console.Write(solution.get(startX + deltaX * j, startY + deltaY * i).ToString("F4") + " ");
                                    }
                                    catch
                                    {
                                        Console.Write("-------- ");
                                    }
                                Console.WriteLine();
                            }

                            Console.WriteLine("Press any key to get back");
                            Console.ReadKey();
                            break;

                        }
                    case 2:
                        {
                            Console.WriteLine("(" + string.Join(" ", solution.coeffs.Select(value => value.ToString("E3"))) + ")");

                            Console.WriteLine("Press any key to get back");
                            Console.ReadKey();
                            break;
                        }
                }
            }
        }
        static Func<double, double, double> ParseExpression(string expr)
        {
            expr = expr.ToLower().Replace("sin", "Math.Sin").Replace("cos", "Math.Cos")
                                 .Replace("tg", "Math.Tan").Replace("ctg", "1.0 / Math.Tan")
                                 .Replace("pow", "Math.Pow").Replace("exp", "Math.Exp")
                                 .Replace("arcsin", "Math.Asin").Replace("arccos", "Math.Acos")
                                 .Replace("arctg", "Math.Atan").Replace("sqrt", "Math.Sqrt")
                                 .Replace("e", "Math.E").Replace("pi", "Math.PI");
            CSharpCodeProvider codeProvider = new CSharpCodeProvider(new Dictionary<string, string>() { { "CompilerVersion", "v3.5" } });
            CompilerParameters parameters = new CompilerParameters(new[] { "System.dll" });
            parameters.GenerateInMemory = true;
            parameters.GenerateExecutable = false;

            CompilerResults results = codeProvider.CompileAssemblyFromSource(parameters,
            @"
            using System;
            class Program 
            {
                public static double func(double x, double y)
                {
                    return " + expr + @";
                }
            }");
            return (Func<double, double, double>)results.CompiledAssembly.GetType("Program").GetMethod("func").CreateDelegate(typeof(Func<double, double, double>));
        }
    }
}