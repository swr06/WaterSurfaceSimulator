#version 450 core 

layout (location = 0) out vec4 o_Color;
layout (location = 1) out vec4 o_Normal;

in vec3 v_WorldPos;

uniform sampler2D u_Texture;

uniform float u_Range;
uniform int u_Res;

layout (std430, binding = 0) buffer SSBO_HM {
	float Heightmap[];
};

int To1DIdx(int x, int y) {
	return (y * u_Res) + x;
}

float Sample(vec2 UV) {
	UV = clamp(UV, 0., 1.);
	ivec2 Texel = ivec2(UV * vec2(float(u_Res)));
	return Heightmap[To1DIdx(Texel.x, Texel.y)];
}

float SamplePx(ivec2 px) {
	return Heightmap[To1DIdx(px.x, px.y)];
}

float Bilinear(vec2 SampleUV)
{
    const ivec2 Cross[4] = ivec2[4](ivec2(0, 0), ivec2(1, 0), ivec2(0, 1), ivec2(1, 1));

    // Relative to center
    vec2 SamplingFragment = (SampleUV * vec2(u_Res)) - vec2(0.5f); 

    // Find how much you need to interpolate across both axis
    vec2 f = fract(SamplingFragment);

    // Fetch 4 neighbours
    float Fetch[4];
    Fetch[0] = SamplePx(ivec2(SamplingFragment) + Cross[0]);
    Fetch[1] = SamplePx(ivec2(SamplingFragment) + Cross[1]);
    Fetch[2] = SamplePx(ivec2(SamplingFragment) + Cross[2]);
    Fetch[3] = SamplePx(ivec2(SamplingFragment) + Cross[3]);

    // Interpolate first based on x position
    float Interp1 = mix(Fetch[0], Fetch[1], float(f.x));
    float Interp2 = mix(Fetch[2], Fetch[3], float(f.x));

    return mix(Interp1, Interp2, float(f.y));
}

float GetHeightAt(vec3 W) {
    vec2 UV = fract(((u_Range + W.xz) / u_Range) * 0.5f);
    return Sample(UV);
}

float GetHeightAtSmooth(vec3 W) {
    vec2 UV = fract(((u_Range + W.xz) / u_Range) * 0.5f);
    return Bilinear(UV);
}

vec3 SampleNormals(vec3 WorldPosition) {

    float Width = (u_Range / float(u_Res));

    vec3 WestVertex = WorldPosition + vec3(Width, 0.0f, 0.0f);
    WestVertex.y = GetHeightAtSmooth(WestVertex);

    vec3 NorthVertex = WorldPosition + vec3(0.0f, 0.0f, Width);
    NorthVertex.y = GetHeightAtSmooth(NorthVertex);

    return normalize(cross(WorldPosition - NorthVertex, WorldPosition - WestVertex));
}


void main() {
	
	vec2 UV = fract(((u_Range + v_WorldPos.xz) / u_Range) * 0.5f);

    float HeightAt = clamp(Bilinear(UV), 0., 1.);
    //HeightAt = pow(HeightAt, 128.0f);

    // Idea : Find vertex normals via screenspace derivatives 
    // Create an orthonormal basis? (TBN), orient High freequency normals to it! 
    vec3 DeltaX = vec3(dFdxFine(v_WorldPos.x), dFdxFine(v_WorldPos.y), dFdxFine(v_WorldPos.z));    
    vec3 DeltaY = vec3(dFdyFine(v_WorldPos.x), dFdyFine(v_WorldPos.y), dFdyFine(v_WorldPos.z));
    vec3 VertexNormal = normalize(cross(DeltaX,DeltaY)); 

    vec2 N = vec2(dFdxFine(HeightAt), dFdyFine(HeightAt)) * 0.5f + 0.5f;

    vec3 Color = vec3(1.);// vec3(0, 117, 119) / 255.;

    Color *= mix(0.1f, 1.0f, clamp(v_WorldPos.y / 3., 0., 1.));

	o_Color = vec4(SampleNormals(vec3(v_WorldPos.x, HeightAt, v_WorldPos.z)),1.0f);

}