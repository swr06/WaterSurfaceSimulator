#version 450 core 

layout (location = 0) out vec4 o_Color;
layout (location = 1) out vec4 o_Normal;

in vec3 v_WorldPos;

uniform sampler2D u_Texture;

uniform float u_Range;
uniform int u_Res;

uniform vec3 u_SunDirection;

uniform bool u_HQ;

layout (std430, binding = 0) buffer SSBO_HM {
	float Heightmap[];
};


float CatmullRom(in vec2 uv);


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

vec2 SpiralPoint(float angle, float scale) {
    return vec2(sin(angle), cos(angle)) * pow(angle / scale, 1.0 / (sqrt(5.0) * 0.5 + 0.5));
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

    if (u_HQ) {
        return CatmullRom(UV);
    }

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

    vec3 Normals = SampleNormals(vec3(v_WorldPos.x, HeightAt, v_WorldPos.z));
    
    o_Color = vec4(Normals,1.);

}



//



float CatmullRom(in vec2 uv)
{
    vec2 texSize = vec2(u_Res);
    vec2 samplePos = uv * texSize;
    vec2 texPos1 = floor(samplePos - 0.5f) + 0.5f;
    vec2 f = samplePos - texPos1;
    vec2 w0 = f * (-0.5f + f * (1.0f - 0.5f * f));
    vec2 w1 = 1.0f + f * f * (-2.5f + 1.5f * f);
    vec2 w2 = f * (0.5f + f * (2.0f - 1.5f * f));
    vec2 w3 = f * f * (-0.5f + 0.5f * f);
    
    vec2 w12 = w1 + w2;
    vec2 offset12 = w2 / (w1 + w2);

    vec2 texPos0 = texPos1 - 1;
    vec2 texPos3 = texPos1 + 2;
    vec2 texPos12 = texPos1 + offset12;

    texPos0 /= texSize;
    texPos3 /= texSize;
    texPos12 /= texSize;

    float result = float(0.0f);

    result +=  Bilinear(vec2(texPos0.x, texPos0.y)) * w0.x * w0.y;
    result +=  Bilinear(vec2(texPos12.x, texPos0.y)) * w12.x * w0.y;
    result +=  Bilinear(vec2(texPos3.x, texPos0.y)) * w3.x * w0.y;
    result +=  Bilinear(vec2(texPos0.x, texPos12.y)) * w0.x * w12.y;
    result +=  Bilinear(vec2(texPos12.x, texPos12.y)) * w12.x * w12.y;
    result +=  Bilinear(vec2(texPos3.x, texPos12.y)) * w3.x * w12.y;
    result +=  Bilinear(vec2(texPos0.x, texPos3.y)) * w0.x * w3.y;
    result +=  Bilinear(vec2(texPos12.x, texPos3.y)) * w12.x * w3.y;
    result +=  Bilinear(vec2(texPos3.x, texPos3.y)) * w3.x * w3.y;

    return result;
}