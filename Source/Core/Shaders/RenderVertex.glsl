#version 330 core

layout (location = 0) in vec2 a_Position;

out vec3 v_WorldPos;

uniform mat4 u_ViewProj;
uniform mat4 u_ViewProjRot;
uniform mat4 u_ModelMatrix;

uniform int u_ResV;
uniform float u_RangeV;

layout (std430, binding = 0) buffer SSBO_HM {
	float Heightmap[];
};

int To1DIdx(int x, int y) {
	return (y * u_ResV) + x;
}

float Sample(vec2 UV) {
	ivec2 Texel = ivec2(fract(UV) * vec2(float(u_ResV)));
	return Heightmap[To1DIdx(Texel.x, Texel.y)];
}

float Sample(ivec2 px) {
	return Heightmap[To1DIdx(px.x, px.y)];
}

float Bilinear(vec2 SampleUV)
{
    const ivec2 Cross[4] = ivec2[4](ivec2(0, 0), ivec2(1, 0), ivec2(0, 1), ivec2(1, 1));

    // Relative to center
    vec2 SamplingFragment = (SampleUV * vec2(u_ResV)) - vec2(0.5f); 

    // Find how much you need to interpolate across both axis
    vec2 f = fract(SamplingFragment);

    // Fetch 4 neighbours
    float Fetch[4];
    Fetch[0] = Sample(ivec2(SamplingFragment) + Cross[0]);
    Fetch[1] = Sample(ivec2(SamplingFragment) + Cross[1]);
    Fetch[2] = Sample(ivec2(SamplingFragment) + Cross[2]);
    Fetch[3] = Sample(ivec2(SamplingFragment) + Cross[3]);

    // Interpolate first based on x position
    float Interp1 = mix(Fetch[0], Fetch[1], float(f.x));
    float Interp2 = mix(Fetch[2], Fetch[3], float(f.x));

    return mix(Interp1, Interp2, float(f.y));
}


void main()
{
	gl_Position = u_ModelMatrix * vec4(a_Position, 1.0f, 1.0f);
	
    v_WorldPos = vec3(gl_Position.xyz);

    vec2 UV = fract(((u_RangeV + v_WorldPos.xz) / u_RangeV) * 0.5f);

    float Height = Bilinear(UV);
    gl_Position.y += Height;
    v_WorldPos.y += Height;

	gl_Position = u_ViewProj * gl_Position;
}