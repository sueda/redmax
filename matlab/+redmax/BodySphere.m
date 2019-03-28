classdef BodySphere < redmax.Body
	%BodySphere A sphere
	%
	
	%%
	properties
		radius
	end
	
	%%
	methods
		%%
		function this = BodySphere(density,radius)
			this = this@redmax.Body(density);
			this.radius = radius;
		end
		
		%%
		function computeInertia_(this)
			% Computes inertia at body
			r = this.radius;
			v = 4/3*pi*r^3;
			m = this.density*v;
			I = 2/5*m*r^2;
			this.I_i = [I I I m m m]';
		end
		
		%%
		function draw_(this)
			n = 8;
			n1 = n+1;
			n1n1 = n1*n1;
			[X,Y,Z] = sphere(n);
			r = this.radius;
			X = r*reshape(X,1,n1n1);
			Y = r*reshape(Y,1,n1n1);
			Z = r*reshape(Z,1,n1n1);
			XYZ = this.E_wi*[X;Y;Z;ones(1,n1n1)];
			X = reshape(XYZ(1,:),n1,n1);
			Y = reshape(XYZ(2,:),n1,n1);
			Z = reshape(XYZ(3,:),n1,n1);
			surf(X,Y,Z,'FaceColor',this.color);
		end
		
		%%
		function s = getAxisSize(this)
			s = 2*this.radius;
		end
	end
end
