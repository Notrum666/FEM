#version 440 core

in vec3 L;
in vec3 coeffs1;
in vec3 coeffs2;

uniform float minTemp;
uniform float maxTemp;

out vec4 outColor;

void main()
{
	float value = coeffs1.x * L[0] * (2f * L[0] - 1f) + coeffs1.y * L[1] * (2f * L[1] - 1f) + coeffs1.z * L[2] * (2f * L[2] - 1f) +
				  4f * (coeffs2.x * L[0] * L[1] + coeffs2.y * L[1] * L[2] + coeffs2.z * L[0] * L[2]);
	value = max(0f, min(1f, (value - minTemp) / (maxTemp - minTemp)));
	outColor = value * vec4(1f, 0f, 0f, 1f) + (1f - value) * vec4(0f, 0f, 1f, 1f);
}  