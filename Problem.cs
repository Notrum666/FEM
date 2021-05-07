using System;
using System.Collections.Generic;
using System.Linq;
using NMath;

namespace FEM
{
    class FirstBoundaryCondition
    {
        public Func<double, double, double, double> function { get; private set; }
        public List<int> points { get; private set; }
        public FirstBoundaryCondition(Func<double, double, double, double> function, List<int> points)
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
        public Func<double, double, double, double> function { get; private set; }
        public List<int> points { get; private set; }
        public SecondBoundaryCondition(Func<double, double, double, double> function, List<int> points)
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
    class Solution
    {
        public Mesh mesh;
        public double[] q, prev_q, prev_prev_q;
        public double from, to, middle;
        public List<KeyValuePair<int, int>>[] adjacencyList;
        public double[] getCoeffs(double t)
        {
            if (t < from || t > to)
                throw new ArgumentOutOfRangeException("t", "Time should be in [from, to] range of current solution.");

            double deltaTime = (from - to) / 2.0;
            double invDoubleDeltaTime = 0.5 / (deltaTime * deltaTime);
            double invDeltaTime = invDoubleDeltaTime * 2.0;
            double[] coeffs = new double[q.Length];
            for (int i = 0; i < q.Length; i++)
                coeffs[i] = q[i] * invDoubleDeltaTime * (t - from) * (t - from + deltaTime) -
                            prev_q[i] * invDeltaTime * (t - from) * (t - to) +
                            prev_prev_q[i] * invDoubleDeltaTime * (t - from + deltaTime) * (t - to);
            return coeffs;
        }
    }
    class Problem
    {
        public Mesh mesh;
        public Func<double, double, double, double> f;
        public Func<double, double, double, double> initialCondition;
        public List<FirstBoundaryCondition> firstBoundaryConditions = new List<FirstBoundaryCondition>();
        public List<SecondBoundaryCondition> secondBoundaryConditions = new List<SecondBoundaryCondition>();
        public double[] lambda;
        public double[] gamma;
        public double deltaTime;

        public Solution solution { get; private set; }

        private int frame = -1;
        private List<KeyValuePair<int, int>>[] adjacencyList;
        private MatrixSymm_Sparse globalMatrix;
        private Vector globalVector;
        private Vector coeffs;
        public void next()
        {
            frame++;
            int size;
            int[] ig;
            int[] jg;
            double[] di;
            double[] gg;
            double[] globalVec;

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
            if (frame == 0)
            {
                int additionalVertexesCount = mesh.vertexes.Length + mesh.polygons.Length - 1;

                size = mesh.vertexes.Length + additionalVertexesCount;
                ig = new int[size + 1];
                jg = new int[mesh.polygons.Length * 6 + additionalVertexesCount * 3];
                di = new double[size];
                gg = new double[mesh.polygons.Length * 6 + additionalVertexesCount * 3];
                globalMatrix = new MatrixSymm_Sparse(size, ig, jg, di, gg);
                globalVector = new Vector(size);
                globalVec = globalVector.values;

                solution = new Solution();
                solution.from = 0;
                solution.middle = 0;
                solution.to = 0;
                solution.mesh = mesh;
                solution.q = new double[size];
                solution.prev_q = new double[size];
                solution.prev_prev_q = new double[size];

                adjacencyList = new List<KeyValuePair<int, int>>[mesh.vertexes.Length];
                for (int i = 0; i < mesh.vertexes.Length; i++)
                    adjacencyList[i] = new List<KeyValuePair<int, int>>();
                solution.adjacencyList = adjacencyList;

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
                                    solution.q[additionalVertexIndex] = initialCondition((mesh.vertexes[i].X + mesh.vertexes[j].X) / 2.0,
                                                                                         (mesh.vertexes[i].Y + mesh.vertexes[j].Y) / 2.0,
                                                                                         0.0);
                                    additionalVertexIndex++;
                                    jg[ig[i + 1]] = j;
                                    ig[i + 1]++;

                                    k = mesh.polygons.Length;
                                    break;
                                }
                        }
                    }
                    solution.q[i] = initialCondition(mesh.vertexes[i].X, mesh.vertexes[i].Y, 0.0);
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

                return;
            }
            else
            {
                size = globalMatrix.size;
                ig = globalMatrix.ig;
                jg = globalMatrix.jg;
                di = globalMatrix.di;
                gg = globalMatrix.gg;
                globalVec = globalVector.values;
                for (int i = 0; i < size; i++)
                {
                    di[i] = 0;
                    globalVec[i] = 0;
                }
                for (int i = 0; i < gg.Length; i++)
                    gg[i] = 0;

                solution.from = solution.middle;
                solution.middle = solution.to;
                solution.to += deltaTime;

                coeffs = new Vector(solution.prev_prev_q);
                solution.prev_prev_q = solution.prev_q;
                solution.prev_q = solution.q;
                solution.q = coeffs.values;
            }

            double[][] localMatrix = new double[6][];
            for (int i = 0; i < 6; i++)
                localMatrix[i] = new double[i + 1];
            double[] localVector = new double[6];
            double[] fValues = new double[6];
            double[] conditionValues = new double[6];
            double[,] alphaMatrix = new double[3, 3];
            double detD;
            Vector2[] v = new Vector2[3];

            double prev_prev_coef = (solution.to - solution.middle) / ((solution.to - solution.from) * (solution.middle - solution.from));
            double prev_coef = (solution.to - solution.from) / ((solution.to - solution.middle) * (solution.middle - solution.from));
            double matrixCoef = (solution.to - solution.from + solution.to - solution.middle) / ((solution.to - solution.from) * (solution.to - solution.middle));
            for (int p = 0; p < mesh.polygons.Length; p++)
            {
                int[] globalIndexes = new int[6];
                for (int i = 0; i < 3; i++)
                {
                    globalIndexes[i] = mesh.polygons[p][i];
                    globalIndexes[i + 3] = getAdditionalVertexBetween(mesh.polygons[p][i], mesh.polygons[p][i == 2 ? 0 : i + 1]);
                }

                v[0] = mesh.vertexes[mesh.polygons[p][0]];
                v[1] = mesh.vertexes[mesh.polygons[p][1]];
                v[2] = mesh.vertexes[mesh.polygons[p][2]];
                detD = (v[1].X - v[0].X) * (v[2].Y - v[0].Y) - (v[2].X - v[0].X) * (v[1].Y - v[0].Y);
                alphaMatrix[0, 0] = (v[1].X * v[2].Y - v[2].X * v[1].Y) / detD;
                alphaMatrix[1, 0] = -(v[0].X * v[2].Y - v[2].X * v[0].Y) / detD;
                alphaMatrix[2, 0] = (v[0].X * v[1].Y - v[1].X * v[0].Y) / detD;
                alphaMatrix[0, 1] = -(v[2].Y - v[1].Y) / detD;
                alphaMatrix[1, 1] = (v[2].Y - v[0].Y) / detD;
                alphaMatrix[2, 1] = -(v[1].Y - v[0].Y) / detD;
                alphaMatrix[0, 2] = (v[2].X - v[1].X) / detD;
                alphaMatrix[1, 2] = -(v[2].X - v[0].X) / detD;
                alphaMatrix[2, 2] = (v[1].X - v[0].X) / detD;

                detD = Math.Abs(detD);

                // mass matrix
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j <= i; j++)
                    {
                        localMatrix[i][j] = (i == j ? 1.0 / 60.0 : -1.0 / 360.0) * detD;
                        localMatrix[i + 3][j + 3] = (i == j ? 4.0 / 45.0 : 2.0 / 45.0) * detD;
                    }
                    localMatrix[3][i] = (i != 2 ? 0 : -1.0 / 90.0) * detD;
                    localMatrix[4][i] = (i != 0 ? 0 : -1.0 / 90.0) * detD;
                    localMatrix[5][i] = (i != 1 ? 0 : -1.0 / 90.0) * detD;
                }

                fValues[0] = f(v[0].X, v[0].Y, solution.to);
                fValues[1] = f(v[1].X, v[1].Y, solution.to);
                fValues[2] = f(v[2].X, v[2].Y, solution.to);
                fValues[3] = f((v[0].X + v[1].X) / 2.0, (v[0].Y + v[1].Y) / 2.0, solution.to);
                fValues[4] = f((v[1].X + v[2].X) / 2.0, (v[1].Y + v[2].Y) / 2.0, solution.to);
                fValues[5] = f((v[0].X + v[2].X) / 2.0, (v[0].Y + v[2].Y) / 2.0, solution.to);
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
                if (frame == 1)
                {
                    double invDeltaTime = 1.0 / deltaTime;
                    for (int i = 0; i < 6; i++)
                    {
                        localMatrix[i][i] *= invDeltaTime;

                        localVector[i] += localMatrix[i][i] * solution.prev_q[globalIndexes[i]];
                        for (int j = 0; j < i; j++)
                        {
                            localMatrix[i][j] *= invDeltaTime;

                            localVector[i] += localMatrix[i][j] * solution.prev_q[globalIndexes[j]];
                            localVector[j] += localMatrix[i][j] * solution.prev_q[globalIndexes[i]];
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < 6; i++)
                    {
                        localVector[i] += localMatrix[i][i] * (prev_coef * solution.prev_q[globalIndexes[i]] - prev_prev_coef * solution.prev_prev_q[globalIndexes[i]]);
                        for (int j = 0; j < i; j++)
                        {
                            localVector[i] += localMatrix[i][j] * (prev_coef * solution.prev_q[globalIndexes[j]] - prev_prev_coef * solution.prev_prev_q[globalIndexes[j]]);
                            localVector[j] += localMatrix[i][j] * (prev_coef * solution.prev_q[globalIndexes[i]] - prev_prev_coef * solution.prev_prev_q[globalIndexes[i]]);
                
                            localMatrix[i][j] *= matrixCoef;
                        }
                        localMatrix[i][i] *= matrixCoef;
                    }
                }

                foreach (SecondBoundaryCondition condition in secondBoundaryConditions)
                {
                    conditionValues[0] = condition.function(v[0].X, v[0].Y, solution.to);
                    conditionValues[1] = condition.function(v[1].X, v[1].Y, solution.to);
                    conditionValues[2] = condition.function(v[2].X, v[2].Y, solution.to);
                    conditionValues[3] = condition.function((v[0].X + v[1].X) / 2.0, (v[0].Y + v[1].Y) / 2.0, solution.to);
                    conditionValues[4] = condition.function((v[1].X + v[2].X) / 2.0, (v[1].Y + v[2].Y) / 2.0, solution.to);
                    conditionValues[5] = condition.function((v[0].X + v[2].X) / 2.0, (v[0].Y + v[2].Y) / 2.0, solution.to);
                    int[] polygon = mesh.polygons[p];
                    int nextPointIndex;
                    double h;
                    for (int i = 0; i < 3; i++)
                    {
                        nextPointIndex = i == 2 ? 0 : i + 1;
                        if (condition.containsEdge(polygon[i], polygon[nextPointIndex]))
                        {
                            h = Math.Sqrt(Math.Pow(v[i].X - v[nextPointIndex].X, 2.0) +
                                          Math.Pow(v[i].Y - v[nextPointIndex].Y, 2.0));
                            localVector[i] += (4.0 * conditionValues[i] +
                                               2.0 * conditionValues[i + 3] -
                                               conditionValues[nextPointIndex]) / 30.0 * h;
                            localVector[nextPointIndex] += (4.0 * conditionValues[nextPointIndex] +
                                               2.0 * conditionValues[i + 3] -
                                               conditionValues[i]) / 30.0 * h;
                            localVector[i + 3] += (2.0 * conditionValues[i] +
                                               16.0 * conditionValues[i + 3] +
                                               2.0 * conditionValues[nextPointIndex]) / 30.0 * h;
                        }
                    }
                }

                // stiffness matrix
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j <= i; j++)
                        localMatrix[i][j] += lambda[p] * (i == j ? 1.0 / 2.0 : -1.0 / 6.0) * (alphaMatrix[i, 1] * alphaMatrix[j, 1] + alphaMatrix[i, 2] * alphaMatrix[j, 2]) * detD;
                    localMatrix[3][i] += lambda[p] * (((i == 1 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[0, 1] + alphaMatrix[i, 2] * alphaMatrix[0, 2])) +
                                                      ((i == 0 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[1, 1] + alphaMatrix[i, 2] * alphaMatrix[1, 2]))) * detD;
                    localMatrix[4][i] += lambda[p] * (((i == 2 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[1, 1] + alphaMatrix[i, 2] * alphaMatrix[1, 2])) +
                                                      ((i == 1 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[2, 1] + alphaMatrix[i, 2] * alphaMatrix[2, 2]))) * detD;
                    localMatrix[5][i] += lambda[p] * (((i == 2 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[0, 1] + alphaMatrix[i, 2] * alphaMatrix[0, 2])) +
                                                      ((i == 0 ? 4.0 / 6.0 : 0) * (alphaMatrix[i, 1] * alphaMatrix[2, 1] + alphaMatrix[i, 2] * alphaMatrix[2, 2]))) * detD;
                }
                localMatrix[3][3] += 4.0 / 3.0 * lambda[p] * (alphaMatrix[0, 1] * alphaMatrix[0, 1] + alphaMatrix[0, 2] * alphaMatrix[0, 2] +
                                                              alphaMatrix[0, 1] * alphaMatrix[1, 1] + alphaMatrix[0, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[1, 1] + alphaMatrix[1, 2] * alphaMatrix[1, 2]) * detD;
                localMatrix[4][4] += 4.0 / 3.0 * lambda[p] * (alphaMatrix[1, 1] * alphaMatrix[1, 1] + alphaMatrix[1, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[2, 1] + alphaMatrix[1, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[2, 1] * alphaMatrix[2, 1] + alphaMatrix[2, 2] * alphaMatrix[2, 2]) * detD;
                localMatrix[5][5] += 4.0 / 3.0 * lambda[p] * (alphaMatrix[2, 1] * alphaMatrix[2, 1] + alphaMatrix[2, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[2, 1] * alphaMatrix[0, 1] + alphaMatrix[2, 2] * alphaMatrix[0, 2] +
                                                              alphaMatrix[0, 1] * alphaMatrix[0, 1] + alphaMatrix[0, 2] * alphaMatrix[0, 2]) * detD;
                localMatrix[4][3] += 2.0 / 3.0 * lambda[p] * (alphaMatrix[0, 1] * alphaMatrix[1, 1] + alphaMatrix[0, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[2, 1] + alphaMatrix[1, 2] * alphaMatrix[2, 2] +
                                                              2.0 * alphaMatrix[0, 1] * alphaMatrix[2, 1] + 2.0 * alphaMatrix[0, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[1, 1] + alphaMatrix[1, 2] * alphaMatrix[1, 2]) * detD;
                localMatrix[5][3] += 2.0 / 3.0 * lambda[p] * (alphaMatrix[0, 1] * alphaMatrix[1, 1] + alphaMatrix[0, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[0, 1] * alphaMatrix[2, 1] + alphaMatrix[0, 2] * alphaMatrix[2, 2] +
                                                              2.0 * alphaMatrix[1, 1] * alphaMatrix[2, 1] + 2.0 * alphaMatrix[1, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[0, 1] * alphaMatrix[0, 1] + alphaMatrix[0, 2] * alphaMatrix[0, 2]) * detD;
                localMatrix[5][4] += 2.0 / 3.0 * lambda[p] * (alphaMatrix[0, 1] * alphaMatrix[2, 1] + alphaMatrix[0, 2] * alphaMatrix[2, 2] +
                                                              alphaMatrix[1, 1] * alphaMatrix[2, 1] + alphaMatrix[1, 2] * alphaMatrix[2, 2] +
                                                              2.0 * alphaMatrix[0, 1] * alphaMatrix[1, 1] + 2.0 * alphaMatrix[0, 2] * alphaMatrix[1, 2] +
                                                              alphaMatrix[2, 1] * alphaMatrix[2, 1] + alphaMatrix[2, 2] * alphaMatrix[2, 2]) * detD;

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
                    setFirstCondition(pointIndex, condition.function(mesh.vertexes[pointIndex].X, mesh.vertexes[pointIndex].Y, solution.to));
                    setFirstCondition(getAdditionalVertexBetween(pointIndex, nextPointIndex),
                        condition.function((mesh.vertexes[pointIndex].X + mesh.vertexes[nextPointIndex].X) / 2.0,
                                           (mesh.vertexes[pointIndex].Y + mesh.vertexes[nextPointIndex].Y) / 2.0,
                                           solution.to));
                }
                setFirstCondition(condition.points[i], condition.function(mesh.vertexes[condition.points[i]].X, mesh.vertexes[condition.points[i]].Y, solution.to));
            }

            SLAE.Solve_CGM_LLt(globalMatrix, globalVector, coeffs, out _);
        }
    }
}
