classdef BodyCylinder < redmax.Body
	%BodyCylinder A cylinder along the Z axis
	%
	
	%%
	properties
		radius
		height
	end
	
	%%
	methods
		%%
		function this = BodyCylinder(density,radius,height)
			this = this@redmax.Body(density);
			this.radius = radius;
			this.height = height;
		end
		
		%%
		function computeInertia_(this)
			% Computes inertia at body
			r2 = this.radius*this.radius;
			h2 = this.height*this.height;
			a = pi*r2;
			v = this.height*a;
			m = v*this.density;
			Ixy = m/12.0*(3*r2 + h2);
			Iz = 0.5*m*r2;
			this.I_i = [Ixy Ixy Iz m m m]';
		end
		
		%%
		function draw_(this)
			n = 8;
			[X,Y,Z] = cylinder(this.radius,n);
			Z = 0.5*(Z*2 - 1)*this.height;
			function [X,Y,Z] = xform(X,Y,Z)
				X = reshape(X,1,2*(n+1));
				Y = reshape(Y,1,2*(n+1));
				Z = reshape(Z,1,2*(n+1));
				XYZ = this.E_wi*[X;Y;Z;ones(1,2*(n+1))];
				X = reshape(XYZ(1,:),2,n+1);
				Y = reshape(XYZ(2,:),2,n+1);
				Z = reshape(XYZ(3,:),2,n+1);
			end
			[X,Y,Z] = xform(X,Y,Z);
			[T,R] = meshgrid(linspace(0,2*pi,n+1),linspace(0,this.radius,2));
			X0 = R.*cos(T);
			Y0 = R.*sin(T);
			Z0 = -0.5*this.height*ones(size(X0));
			[X_,Y_,Z_] = xform(X0,Y0,Z0);
			X = [X,X_];
			Y = [Y,Y_];
			Z = [Z,Z_];
			Z0 = 0.5*this.height*ones(size(X0));
			[X_,Y_,Z_] = xform(X0,Y0,Z0);
			X = [X,X_];
			Y = [Y,Y_];
			Z = [Z,Z_];
			surf(X,Y,Z,'FaceColor',this.color);
		end
		
		%%
		function s = getAxisSize(this)
			s = 2*this.radius;
		end
	end
end
