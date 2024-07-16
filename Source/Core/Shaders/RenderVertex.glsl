#version 450 core

layout (location = 0) in vec2 a_Position;

out vec3 v_WorldPos;

uniform mat4 u_ViewProj;
uniform mat4 u_ViewProjRot;
uniform mat4 u_ModelMatrix;

uniform int u_ResV;
uniform float u_RangeV;

uniform bool u_HQv;

layout (std430, binding = 0) buffer SSBO_HM {
	float Heightmap[];
};

int To1DIdx(int x, int y) {
	return (y * u_ResV) + x;
}

float CatmullRom(in vec2 uv);
float Bicubic(vec2 coord);

float Sample(vec2 UV) {
	ivec2 Texel = ivec2((UV) * vec2(float(u_ResV)));
	return Heightmap[To1DIdx(Texel.x, Texel.y)];
}

float Sample(ivec2 px) {
    
    px = clamp(px, ivec2(0), ivec2(u_ResV - 1));
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

float SampleHeight(vec2 UV){
    bool HQ = u_HQv;
    if (HQ) {
        return CatmullRom(UV);
    }

    return Bilinear(UV);
}


void main()
{
	gl_Position = u_ModelMatrix * vec4(a_Position * u_RangeV, 1.0f, 1.0f);
    gl_Position.y += 1.0f;
	
    v_WorldPos = vec3(gl_Position.xyz);

    vec2 UV = (((u_RangeV + v_WorldPos.xz) / u_RangeV) * 0.5f);

    float Height = SampleHeight(UV);
    gl_Position.y += Height;
    v_WorldPos.y += Height;

	gl_Position = u_ViewProj * gl_Position;
}

float CatmullRom(in vec2 uv)
{
    vec2 texSize = vec2(u_ResV);
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

vec4 cubic(float x) {
  float x2 = x * x;
  float x3 = x2 * x;
  vec4 w;
  w.x =   -x3 + 3*x2 - 3*x + 1;
  w.y =  3*x3 - 6*x2       + 4;
  w.z = -3*x3 + 3*x2 + 3*x + 1;
  w.w =  x3;
  return w / 6.f;
}

float Bicubic(vec2 coord) {
    coord *= float(u_ResV);

    float fx = fract(coord.x);
    float fy = fract(coord.y);
    coord.x -= fx;
    coord.y -= fy;

    vec4 xcubic = cubic(fx);
    vec4 ycubic = cubic(fy);

    vec4 c = vec4(coord.x - 0.5, coord.x + 1.5, coord.y - 0.5, coord.y + 1.5);
    vec4 s = vec4(xcubic.x + xcubic.y, xcubic.z + xcubic.w, ycubic.x + ycubic.y, ycubic.z + ycubic.w);
    vec4 offset = c + vec4(xcubic.y, xcubic.w, ycubic.y, ycubic.w) / s;

    float sample0 = Bilinear(vec2(offset.x, offset.z) / float(u_ResV));
    float sample1 = Bilinear(vec2(offset.y, offset.z) / float(u_ResV));
    float sample2 = Bilinear(vec2(offset.x, offset.w) / float(u_ResV));
    float sample3 = Bilinear(vec2(offset.y, offset.w) / float(u_ResV));

    float sx = s.x / (s.x + s.y);
    float sy = s.z / (s.z + s.w);

    return mix( mix(sample3, sample2, sx), mix(sample1, sample0, sx), sy);
}