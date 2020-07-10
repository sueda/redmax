classdef JointRevolute < redmax.Joint
	%%
	properties
		axis     % Axis of rotation
		radius   % Radius for contact
		height   % Height for contact
	end
	
	%%
	methods
		%%
		function this = JointRevolute(parent,body,axis)
			this = this@redmax.Joint(parent,body,1);
			this.axis = axis/norm(axis);
			this.radius = 1.0;
			this.height = 1.0;
		end
		
		%%
		function setGeometry(this,radius,height)
			this.radius = radius;
			this.height = height;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this)
			n = this.ndof;
			q = this.q;
			qdot = this.qdot;
			a = this.axis;
			
			this.Q = eye(4);
			R = se3.aaToMat(a,q); % faster than se3.exp(se3.brac(a*q));
			this.Q(1:3,1:3) = R;
			this.A = se3.Ad(this.Q);
			this.S = [a; zeros(3,1)];
			
			this.Sdot = zeros(6,n);
			this.dSdotdq = zeros(6,n,n);
			
			this.dSdq = zeros(6,n,n);
			this.dAdq = zeros(6,6,n);
			abrac = se3.brac(a);
			dRdq = R*abrac;
			this.dAdq(1:3,1:3) = dRdq;
			this.dAdq(4:6,4:6) = dRdq;
			
			Rdot = dRdq*qdot;
			this.Adot = zeros(6,6);
			this.Adot(1:3,1:3) = Rdot;
			this.Adot(4:6,4:6) = Rdot;
			d2Rdq2 = dRdq*abrac;
			tmp = d2Rdq2*qdot;
			this.dAdotdq = zeros(6,6,n);
			this.dAdotdq(1:3,1:3) = tmp;
			this.dAdotdq(4:6,4:6) = tmp;
		end
		
		%%
		function draw_(this)
			n = 8;
			E_ja = eye(4);
			z = [0 0 1]';
			angle = acos(max(-1.0, min(this.axis'*z, 1.0)));
			E_ja(1:3,1:3) = se3.aaToMat(cross(this.axis,z),angle);
			E = this.E_wj*E_ja;
			% Cylinder
			r = this.radius;
			h = this.height;
			S = diag([1,1,h,1]);
			T = [eye(3),[0 0 -h/2]'; 0 0 0 1];
			[X,Y,Z] = cylinder(r,n);
			X = reshape(X,1,2*(n+1));
			Y = reshape(Y,1,2*(n+1));
			Z = reshape(Z,1,2*(n+1));
			XYZ = E*T*S*[X;Y;Z;ones(1,2*(n+1))];
			X0 = reshape(XYZ(1,:),2,n+1);
			Y0 = reshape(XYZ(2,:),2,n+1);
			Z0 = reshape(XYZ(3,:),2,n+1);
			%surf(X,Y,Z,'FaceColor',this.body.color);
			% Ends
			% https://www.mathworks.com/matlabcentral/answers/320174-plot-a-surface-on-a-circle-and-annulus
			[T,R] = meshgrid(linspace(0,2*pi,n+1),linspace(0,r,2));
			X = R.*cos(T);
			Y = R.*sin(T);
			Z = -h/2*ones(size(R));
			X = reshape(X,1,2*(n+1));
			Y = reshape(Y,1,2*(n+1));
			Z = reshape(Z,1,2*(n+1));
			XYZ = E*[X;Y;Z;ones(1,2*(n+1))];
			X1 = reshape(XYZ(1,:),2,n+1);
			Y1 = reshape(XYZ(2,:),2,n+1);
			Z1 = reshape(XYZ(3,:),2,n+1);
			%surf(X,Y,Z,'FaceColor',this.body.color);
			X = R.*cos(T);
			Y = R.*sin(T);
			Z = h/2*ones(size(R));
			X = reshape(X,1,2*(n+1));
			Y = reshape(Y,1,2*(n+1));
			Z = reshape(Z,1,2*(n+1));
			XYZ = E*[X;Y;Z;ones(1,2*(n+1))];
			X2 = reshape(XYZ(1,:),2,n+1);
			Y2 = reshape(XYZ(2,:),2,n+1);
			Z2 = reshape(XYZ(3,:),2,n+1);
			surf([X0 X1 X2],[Y0 Y1 Y2],[Z0 Z1 Z2],'FaceColor',this.body.color);
		end
	end
end
