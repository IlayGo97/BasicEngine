#version 330
uniform mat4 MVP;
uniform mat4 Normal;

in vec3 position;
in vec3 color;
in vec3 normal;
in vec2 texCoords;


out vec3 color0;
out vec3 normal0;

void main()
{
	color0 = color;
	normal0 = vec3(Normal*vec4(normal,0));
	gl_Position = MVP *Normal* vec4(position, 1.0);
}
