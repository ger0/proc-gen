#version 330

/* fragment shader for objecs with many meshes with different textures */

uniform sampler2D tex; // texturing unit

out vec4 pixelColor;

in vec4 i_color;
in vec4 n;
in vec4 l[2];
in vec4 v;
in vec2 i_texc;

vec4 calc_light(vec4 light_position, vec4 normal, vec4 vertex, vec4 surface_color, vec4 reflect_color){

	vec4 ml = normalize(light_position);
	vec4 mr = reflect(-ml, normal); //w przestrzeni oka, bo l i n są w tej przestrzeni

	float nl = clamp(dot(normal, ml), 0,1); // cos kąta
	float alfa = 25;
	float rv = pow(clamp(dot(mr,vertex), 0, 1), alfa);

	return vec4(surface_color.rgb * nl, surface_color.a) + vec4(reflect_color.rgb*rv, 0);
}


void main(void) {
	float brighteness = 1.5;

	vec4 kd = texture(tex,i_texc) * brighteness; //kolor powierzchni
	vec4 mn = normalize(n);
	vec4 mv = normalize(v);
	vec4 ks = kd;	//kolor światła odbitego


//L = ka*la + kd*ld*nl + ks*ls*pow(rv, alfa)
	// pixelColor=vec4(kd.rgb * nl, kd.a) + vec4(ks.rgb*rv, 0);

	pixelColor = vec4(0,0,0,0);
	for (int i=0; i<2; i++){
		pixelColor += calc_light(l[i], mn, mv, kd, ks);
	}



}
