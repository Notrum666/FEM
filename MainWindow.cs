using System;
using System.Collections.Generic;
using System.Linq;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL4;
using System.IO;
using Microsoft.CSharp;
using System.CodeDom.Compiler;
using OpenTK.Input;
using System.Globalization;

namespace FEM
{
    public class MainWindow : GameWindow
    {
        private Problem problem;
        private double time = 0;
        private Camera camera;
        private Shader shader;
        private double minTemp, maxTemp;
        private Func<double, double, double, double> deltaTime;
        private bool disableLimit;
        private double timeScale;
        public MainWindow(int width, int height) : base(width, height, GraphicsMode.Default, "FEM", GameWindowFlags.FixedWindow)
        {
            Input.Init();

            GL.Disable(EnableCap.DepthTest);
            GL.Enable(EnableCap.CullFace);
            GL.CullFace(CullFaceMode.Back);

            camera = new Camera(16f, 9f);

        }
        protected override void OnLoad(EventArgs e)
        {
            problem = new Problem();
            problem.mesh = Mesh.Parse("mesh.txt");
            string[] materials = File.ReadAllLines("materials.txt");
            problem.lambda = materials[0].Split(' ').Select(value => double.Parse(value)).ToArray();
            problem.gamma = materials[1].Split(' ').Select(value => double.Parse(value)).ToArray();
            
            problem.f = ParseExpression(File.ReadAllText("function.txt"));
            string[] lines = File.ReadAllLines("first boundary.txt");
            for (int i = 0; i < lines.Length; i += 2)
                problem.firstBoundaryConditions.Add(new FirstBoundaryCondition(ParseExpression(lines[i]), lines[i + 1].Split(' ').Select(value => int.Parse(value)).ToList()));
            lines = File.ReadAllLines("second boundary.txt");
            for (int i = 0; i < lines.Length; i += 2)
                problem.secondBoundaryConditions.Add(new SecondBoundaryCondition(ParseExpression(lines[i]), lines[i + 1].Split(' ').Select(value => int.Parse(value)).ToList()));
            lines = File.ReadAllLines("initial condition.txt");
            problem.initialCondition = ParseExpression(lines[0]);
            lines = File.ReadAllLines("config.txt");
            deltaTime = (x, y, t) => 1.0;
            minTemp = 0.0;
            maxTemp = 1.0;
            timeScale = 1.0;
            disableLimit = false;
            foreach (string line in lines)
            {
                string[] words = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                switch (words[0].ToLower())
                {
                    case "dt":
                    case "deltatime":
                        deltaTime = ParseExpression(words[1]); //double.Parse(words[1].Replace(',', '.'), CultureInfo.InvariantCulture);
                        break;
                    case "mintemp":
                        minTemp = double.Parse(words[1].Replace(',', '.'), CultureInfo.InvariantCulture);
                        break;
                    case "maxtemp":
                        maxTemp = double.Parse(words[1].Replace(',', '.'), CultureInfo.InvariantCulture);
                        break;
                    case "disablelimit":
                    case "dislimit":
                        disableLimit = true;
                        break;
                    case "timescale":
                    case "scale":
                        timeScale = double.Parse(words[1].Replace(',', '.'), CultureInfo.InvariantCulture);
                        break;
                }
            }

            problem.deltaTime = deltaTime(0, 0, 0);
            problem.next();
            problem.deltaTime = deltaTime(0, 0, problem.solution.to);
            problem.next();
            problem.deltaTime = deltaTime(0, 0, problem.solution.to);
            problem.next();

            shader = AssetsManager.LoadShader("default", new ShaderComponent("Assets\\Shaders\\default.fsh"),
                                                         new ShaderComponent("Assets\\Shaders\\default.gsh"),
                                                         new ShaderComponent("Assets\\Shaders\\default.vsh"));
        }
        int getAdditionalVertexBetween(int i, int j)
        {
            if (j > i)
            {
                i = i + j;
                j = i - j;
                i = i - j;
            }
            foreach (KeyValuePair<int, int> pair in problem.solution.adjacencyList[i])
                if (pair.Key == j)
                    return pair.Value;
            return -1;
        }
        protected override void OnRenderFrame(FrameEventArgs e)
        {
            time += e.Time * timeScale;
            while (time > problem.solution.to)
            {
                problem.deltaTime = deltaTime(0, 0, problem.solution.to);
                if (!disableLimit)
                    problem.deltaTime = Math.Max(0.01, problem.deltaTime);
                problem.next();
            }

            GL.Clear(ClearBufferMask.ColorBufferBit);
            GL.ClearColor(0.2f, 0.2f, 0.2f, 1.0f);

            GL.UseProgram(shader.id);

            double[] coeffs = problem.solution.getCoeffs(time);
            float[] data = new float[problem.mesh.polygons.Length * 12];
            int curArrayIndex = 0;
            foreach (int[] poly in problem.mesh.polygons)
            {
                for (int j = 0; j < 3; j++, curArrayIndex += 4)
                {
                    data[curArrayIndex] = (float)problem.mesh.vertexes[poly[j]].X;
                    data[curArrayIndex + 1] = (float)problem.mesh.vertexes[poly[j]].Y;
                    data[curArrayIndex + 2] = (float)coeffs[poly[j]];
                    data[curArrayIndex + 3] = (float)coeffs[getAdditionalVertexBetween(poly[j], j != 2 ? poly[j + 1] : poly[0])];
                }
            }

            Matrix4 camSpace = camera.camSpace;
            GL.UniformMatrix4(shader.locations["camSpace"], true, ref camSpace);
            GL.Uniform1(shader.locations["maxTemp"], (float)maxTemp);
            GL.Uniform1(shader.locations["minTemp"], (float)minTemp);

            GL.BindVertexArray(problem.mesh.VAO);
            GL.BindBuffer(BufferTarget.ArrayBuffer, problem.mesh.VBO);
            GL.BufferData(BufferTarget.ArrayBuffer, sizeof(float) * data.Length, data, BufferUsageHint.StaticDraw);
            GL.DrawArrays(PrimitiveType.Triangles, 0, problem.mesh.polygons.Length * 3);
            GL.BindBuffer(BufferTarget.ArrayBuffer, 0);
            GL.BindVertexArray(0);

            SwapBuffers();
        }
        protected override void OnUpdateFrame(FrameEventArgs e)
        {
            Input.OnUpdateFrame();

            float scroll = Input.GetScrollDelta();
            if (scroll != 0)
                camera.scale *= scroll > 0 ? 1.1f : 0.9f;
            double cameraSpeed = 10.0;
            if (Input.IsKeyDown(Key.Right))
                camera.position.X += (float)(cameraSpeed * e.Time) / camera.scale;
            if (Input.IsKeyDown(Key.Left))
                camera.position.X -= (float)(cameraSpeed * e.Time) / camera.scale;
            if (Input.IsKeyDown(Key.Up))
                camera.position.Y += (float)(cameraSpeed * e.Time) / camera.scale;
            if (Input.IsKeyDown(Key.Down))
                camera.position.Y -= (float)(cameraSpeed * e.Time) / camera.scale;
            if (Input.IsKeyPressed(Key.Escape))
                Exit();
        }
        static Func<double, double, double, double> ParseExpression(string expr)
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
                public static double func(double x, double y, double t)
                {
                    return " + expr + @";
                }
            }");
            foreach (CompilerError error in results.Errors)
                if (!error.IsWarning)
                    throw new Exception(expr + " is in invalid format.");
            return (Func<double, double, double, double>)results.CompiledAssembly.GetType("Program").GetMethod("func").CreateDelegate(typeof(Func<double, double, double, double>));
        }
    }
}
