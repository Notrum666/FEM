#version 440 core
layout(location = 0) in vec2 v;
layout(location = 1) in vec2 coeffs;

uniform mat4 camSpace;

out vec2 _coeffs;

void main()
{
    gl_Position = camSpace * vec4(v, 0.0f, 1.0f);
    _coeffs = coeffs;
}