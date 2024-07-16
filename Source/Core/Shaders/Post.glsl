#version 330 core 

#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))


layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_Texture;

uniform float u_Exposure;

uniform float u_Time;
uniform bool u_Distort;

float perlin_noise(vec2 uv, float cells_count);
vec3 blur(sampler2D sp, vec2 uv, vec2 scale);

float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}


// ACES Tonemap operator 
mat3 ACESInputMat = mat3(
    0.59719, 0.07600, 0.02840,
    0.35458, 0.90834, 0.13383,
    0.04823, 0.01566, 0.83777
);

// ODT_SAT => XYZ => D60_2_D65 => sRGB
mat3 ACESOutputMat = mat3(
    1.60475, -0.10208, -0.00327,
    -0.53108, 1.10813, -0.07276,
    -0.07367, -0.00605, 1.07602
);

vec3 RRTAndODTFit(vec3 v)
{
    vec3 a = v * (v + 0.0245786f) - 0.000090537f;
    vec3 b = v * (0.983729f * v + 0.4329510f) + 0.238081f;
    return a / b;
}

vec3 ACESFitted(vec3 Color)
{ 
    Color.rgb = ACESInputMat * Color.rgb;
    Color.rgb = RRTAndODTFit(Color.rgb);
    Color.rgb = ACESOutputMat * Color.rgb;

    return clamp(Color, 0., 1.);
}

vec3 ToSRGB(vec3 x) {
    return mix(1.055*pow(x, vec3(1./2.4)) - 0.055, 12.92*x, step(x,vec3(0.0031308)));
}



const vec4 u_WaveStrengthX = vec4(4.15, 4.66, 0.0016, 0.0015);
const vec4 u_WaveStrengthY = vec4(2.54, 6.33, 0.00102, 0.0025);

vec2 WaterDistort(vec2 TextureCoordinates) {

    float Hash = 16. + perlin_noise(mod(v_TexCoords*0.4f + vec2(u_Time*1.5f, u_Time*0.8f), 64.0f),10.0f);

    // Distort TextureCoordinates coordinates in the y-direction
    TextureCoordinates.y += (cos((TextureCoordinates.y + u_Time * u_WaveStrengthY.y + u_WaveStrengthY.x * Hash)) * u_WaveStrengthY.z) +
            (cos((TextureCoordinates.y + u_Time) * 10.0) * u_WaveStrengthY.w);

    // Distort TextureCoordinates coordinates in the x-direction
    TextureCoordinates.x += (sin((TextureCoordinates.y + u_Time * u_WaveStrengthX.y + u_WaveStrengthX.x * Hash)) * u_WaveStrengthX.z) +
            (sin((TextureCoordinates.y + u_Time) * 15.0) * u_WaveStrengthX.w);

    if (TextureCoordinates != clamp(TextureCoordinates,0.,1.)) {
        return v_TexCoords;
    }

    return TextureCoordinates;
}

void main() {

    vec2 g_TexCoords = u_Distort?WaterDistort(v_TexCoords):v_TexCoords;
    float Hash = bayer16(gl_FragCoord.xy);
    float Exposure = u_Exposure;
    vec3 Sample = texture(u_Texture, g_TexCoords).xyz;
    vec2 TexelSize = 1./vec2(textureSize(u_Texture,0));
	o_Color.xyz = (ACESFitted(Sample * Exposure));
    o_Color.w = 1.;

    vec2 TexSize = textureSize(u_Texture, 0).xy;

    float Aspect = TexSize.x / TexSize.y;

}





// By spectrespect
vec2 n22 (vec2 p)
{
    vec3 a = fract(p.xyx * vec3(123.34, 234.34, 345.65));
    a += dot(a, a + 34.45);
    return fract(vec2(a.x * a.y, a.y * a.z));
}

vec2 get_gradient(vec2 pos)
{
    float twoPi = 6.283185;
    float angle = n22(pos).x * twoPi;
    return vec2(cos(angle), sin(angle));
}

float perlin_noise(vec2 uv, float cells_count)
{
    vec2 pos_in_grid = uv * cells_count;
    vec2 cell_pos_in_grid =  floor(pos_in_grid);
    vec2 local_pos_in_cell = (pos_in_grid - cell_pos_in_grid);
    vec2 blend = local_pos_in_cell * local_pos_in_cell * (3.0f - 2.0f * local_pos_in_cell);
    
    vec2 left_top = cell_pos_in_grid + vec2(0, 1);
    vec2 right_top = cell_pos_in_grid + vec2(1, 1);
    vec2 left_bottom = cell_pos_in_grid + vec2(0, 0);
    vec2 right_bottom = cell_pos_in_grid + vec2(1, 0);
    
    float left_top_dot = dot(pos_in_grid - left_top, get_gradient(left_top));
    float right_top_dot = dot(pos_in_grid - right_top,  get_gradient(right_top));
    float left_bottom_dot = dot(pos_in_grid - left_bottom, get_gradient(left_bottom));
    float right_bottom_dot = dot(pos_in_grid - right_bottom, get_gradient(right_bottom));
    
    float noise_value = mix(
                            mix(left_bottom_dot, right_bottom_dot, blend.x), 
                            mix(left_top_dot, right_top_dot, blend.x), 
                            blend.y);
   
    
    return (0.5 + 0.5 * (noise_value / 0.7));
}

#define pow2(x) (x * x)

const float pi = atan(1.0) * 4.0;
const int samples = 8;
const float sigma = float(samples) * 0.25;

float gaussian(vec2 i) {
    return 1.0 / (2.0 * pi * pow2(sigma)) * exp(-((pow2(i.x) + pow2(i.y)) / (2.0 * pow2(sigma))));
}

vec3 blur(sampler2D sp, vec2 uv, vec2 scale) {
    vec3 col = vec3(0.0);
    float accum = 0.0;
    float weight;
    vec2 offset;
    
    for (int x = -samples / 2; x < samples / 2; ++x) {
        for (int y = -samples / 2; y < samples / 2; ++y) {
            offset = vec2(x, y);
            weight = gaussian(offset);
            col += texture(sp, uv + scale * offset).rgb * weight;
            accum += weight;
        }
    }
    
    return col / accum;
}