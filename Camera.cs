using System;
using OpenTK;

namespace FEM
{
    public class Camera
    {
        public Vector2 position = Vector2.Zero;
        public float width, height, scale = 1.0f;
        public Matrix4 camSpace
        {
            get
            {
                return new Matrix4(scale / width, 0f, 0f, -scale * position.X / width,
                                   0f, scale / height, 0f, -scale * position.Y / height,
                                   0f, 0f, 1f, -0.5f,
                                   0f, 0f, 0f, 1f);
            }
        }
        public Camera(float width, float height)
        {
            this.width = width;
            this.height = height;
        }
    }
}
