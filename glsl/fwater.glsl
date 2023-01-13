#version 330
precision mediump float;

in vec3  iNormal;
in vec3  vertPos;
in vec4  iColor;
in float visibility;

uniform vec3 cameraPos;
uniform vec3 skyColor;
uniform sampler2D depthTexture;

uniform float near;
uniform float far;

uniform vec3 sunDir;
const vec3 lightColor = vec3(0.6, 0.6, 0.5);
const float screenGamma = 2.2; // Assume the monitor is calibrated to the sRGB color space

float ld(float original_depth) {
    return (2.0 * near) / (far + near - original_depth * (far - near));
}

void main() {
	vec3 ambientColor = lightColor * 0.1;
	vec3 normal = normalize(iNormal);

	vec3 cameraDir = cameraPos - vertPos;
	cameraDir = normalize(cameraDir);

	float lambertian = max(dot(sunDir, normal), 0.0);
	float lambertian1 = max(dot(cameraDir, normal), 0.0);

	vec3 diffuse = lambertian * lightColor;

	vec3 viewDir = normalize(cameraPos - vertPos);
	vec3 reflectDir = reflect(-sunDir, normal);  

	float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
	vec3 specular = 0.5 * spec * lightColor;

	vec3 colorLinear = (ambientColor + diffuse + specular) * iColor.rgb;

	vec3 colorCorrect = pow(colorLinear, vec3(1.0 / screenGamma));
	vec3 outColor = mix(skyColor, colorCorrect, visibility);
	//gl_FragColor = vec4(outColor, iColor.a);
	float olddepth = texelFetch(depthTexture, ivec2(gl_FragCoord.xy), 0).r;
	//float alpha = depth;
	float depth = pow(ld(olddepth) - ld(gl_FragDepth), 0.99) + 0.15;
	float alpha = clamp(depth, 0.16, 0.9);
	gl_FragColor = vec4(outColor, alpha);
	//gl_FragColor = vec4(gl_FragCoord.x / 800, 0, 0, 1.0);
}
