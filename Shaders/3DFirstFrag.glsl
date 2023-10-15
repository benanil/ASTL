#version 330

in  vec3 vnorm;
// in  vec2 vtexCoord;

out vec4 fragColor;

void main()
{
    // float ndl = max(dot(vnorm, vec3(sin(vtexCoord.x), cos(vtexCoord.y), 0.0f)), 0.2);
    fragColor = vec4(vnorm,  1.0f);
}