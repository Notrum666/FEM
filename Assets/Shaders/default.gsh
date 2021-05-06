#version 440 core

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

in vec2 _coeffs[];

out vec3 L;
out vec3 coeffs1;
out vec3 coeffs2;

void main()
{
    coeffs1 = vec3(_coeffs[0].x, _coeffs[1].x, _coeffs[2].x);
    coeffs2 = vec3(_coeffs[0].y, _coeffs[1].y, _coeffs[2].y);

    gl_Position = gl_in[0].gl_Position;
    L = vec3(1f, 0f, 0f);
    EmitVertex();

    gl_Position = gl_in[1].gl_Position;
    L = vec3(0f, 1f, 0f);
    EmitVertex();

    gl_Position = gl_in[2].gl_Position;
    L = vec3(0f, 0f, 1f);
    EmitVertex();

    EndPrimitive();
}