using System;
using System.Collections.Generic;
using NMath;
using System.IO;
using OpenTK.Graphics.OpenGL4;

namespace FEM
{
    class Mesh
    {
        public Vector2[] vertexes { get; private set; }
        public int[][] polygons { get; private set; }
        public int VAO { get; private set; }
        public int VBO { get; private set; }
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

            mesh.VAO = GL.GenVertexArray();
            GL.BindVertexArray(mesh.VAO);

            mesh.VBO = GL.GenBuffer();
            GL.BindBuffer(BufferTarget.ArrayBuffer, mesh.VBO);

            GL.VertexAttribPointer(0, 2, VertexAttribPointerType.Float, false, 4 * 4, 0);
            GL.VertexAttribPointer(1, 2, VertexAttribPointerType.Float, false, 4 * 4, 2 * 4);

            GL.EnableVertexAttribArray(0);
            GL.EnableVertexAttribArray(1);

            GL.BindBuffer(BufferTarget.ArrayBuffer, 0);

            GL.BindVertexArray(0);

            return mesh;
        }
    }
}
