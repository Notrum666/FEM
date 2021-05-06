using System;
using System.Collections.Generic;
using System.IO;
using OpenTK.Graphics.OpenGL4;

namespace FEM
{
    public class Shader
    {
        public int id;
        public Dictionary<string, int> locations;
        public Shader(int id)
        {
            this.id = id;
            locations = new Dictionary<string, int>();
            int uniformsCount;
            GL.GetProgram(id, GetProgramParameterName.ActiveUniforms, out uniformsCount);
            for (int i = 0; i < uniformsCount; i++)
            {
                string uniformName = GL.GetActiveUniform(id, i, out _, out _);
                if (uniformName.EndsWith("[0]"))
                    uniformName = uniformName.Substring(0, uniformName.Length - 3);
                locations[uniformName] = GL.GetUniformLocation(id, uniformName);
            }
        }
    }
    public class ShaderComponent
    {
        public int Id { get; private set; }
        public ShaderType Type { get; private set; }
        public ShaderComponent(string path)
        {
            if (!Directory.Exists(Path.GetDirectoryName(path)) || !File.Exists(path))
                throw new FileNotFoundException("Vertex shader file not found", path);

            switch (Path.GetExtension(path))
            {
                case ".vsh":
                    Type = ShaderType.VertexShader;
                    break;
                case ".fsh":
                    Type = ShaderType.FragmentShader;
                    break;
                case ".gsh":
                    Type = ShaderType.GeometryShader;
                    break;
                case ".csh":
                    Type = ShaderType.ComputeShader;
                    break;
                default:
                    throw new ArgumentException("Unable to get shader type from file extension, change the file extension or define shader type explicitly by using other constructor overload.");
            }
            Id = GL.CreateShader(Type);
            GL.ShaderSource(Id, File.ReadAllText(path));
            GL.CompileShader(Id);
            int result;
            GL.GetShader(Id, ShaderParameter.CompileStatus, out result);
            if (result == 0)
                throw new Exception("Shader compilation error, shader type: " + Type.ToString() + ", error: " + GL.GetShaderInfoLog(Id));
        }
        public ShaderComponent(ShaderType type, string path)
        {
            if (!Directory.Exists(Path.GetDirectoryName(path)) || !File.Exists(path))
                throw new FileNotFoundException("Vertex shader file not found", path);

            Type = type;
            Id = GL.CreateShader(type);
            GL.ShaderSource(Id, File.ReadAllText(path));
            GL.CompileShader(Id);
            int result;
            GL.GetShader(Id, ShaderParameter.CompileStatus, out result);
            if (result == 0)
                throw new Exception("Shader compilation error, shader type: " + Type.ToString() + ", error: " + GL.GetShaderInfoLog(Id));
        }
    }
    public static class AssetsManager
    {
        public static Dictionary<string, Shader> Shaders = new Dictionary<string, Shader>();
        public static Shader LoadShader(string shaderName, params ShaderComponent[] shaderComponents)
        {
            int program = GL.CreateProgram();
            foreach (ShaderComponent component in shaderComponents)
                GL.AttachShader(program, component.Id);
            GL.LinkProgram(program);
            int result;
            GL.GetProgram(program, GetProgramParameterName.LinkStatus, out result);
            if (result == 0)
                throw new Exception("Program linking error: " + GL.GetProgramInfoLog(program));

            foreach (ShaderComponent component in shaderComponents)
            {
                GL.DetachShader(program, component.Id);
                GL.DeleteShader(component.Id);
            }

            Shader shader = new Shader(program);
            Shaders[shaderName] = shader;
            return shader;
        }
    }
}
