classdef ConstraintLoop < redmax.Constraint
	% ConstraintLoop Spherical constraint between two bodies
	
	%%
	properties
		bodyA
		bodyB
		xA    % Local position wrt A
		xB    % Local position wrt B
		
		% These variables are for revolute constraints only. They should be
		% moved to a subclass
		radius
		height
	end
	
	%%
	methods
		%%
		function this = ConstraintLoop(bodyA,bodyB)
			this = this@redmax.Constraint(2,0,0,0);
			this.name = [bodyA.name,'-',bodyB.name];
			this.bodyA = bodyA;
			this.bodyB = bodyB;
			this.radius = 1.0;
			this.height = 1.0;
		end
		
		%%
		function setPositions(this,xA,xB)
			this.xA = xA;
			this.xB = xB;
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
		function [Gm,Gmdot,gm,gmdot,gmddot] = computeJacEqM_(this,Gm,Gmdot,gm,gmdot,gmddot)
			rows = this.idxEM;
			E_wa = this.bodyA.E_wi;
			E_wb = this.bodyB.E_wi;
			R_wa = E_wa(1:3,1:3);
			R_wb = E_wb(1:3,1:3);
			% Get two directions orthonormal to A's hinge axis
			jointA = this.bodyA.joint;
			v0 = R_wa*jointA.axis;
			[~,imin] = min(abs(v0));
			v1 = zeros(3,1);
			v1(imin) = 1;
			v2 = cross(v0,v1);
			v2 = v2/norm(v2);
			v1 = cross(v2,v0);
			v1 = v1/norm(v1);
			v12 = [v1 v2];
			GammaA = se3.Gamma(this.xA);
			GammaB = se3.Gamma(this.xB);
			waBrac = se3.brac(this.bodyA.phi(1:3));
			wbBrac = se3.brac(this.bodyB.phi(1:3));
			colsA = this.bodyA.idxM;
			colsB = this.bodyB.idxM;
			this.idxQ = [colsA colsB];
			Gm(rows,colsA) =  v12'*R_wa*GammaA;
			Gm(rows,colsB) = -v12'*R_wb*GammaB;
			Gmdot(rows,colsA) =  v12'*R_wa*waBrac*GammaA;
			Gmdot(rows,colsB) = -v12'*R_wb*wbBrac*GammaB;
			dx = E_wa*[this.xA;1] - E_wb*[this.xB;1];
			gm(rows) = v12'*dx(1:3);
		end
		
		%%
		function generateContactsJoint_(this)
			% Assumes a revolute constraint!
			% The contact points will be generated in the 'C' frame, which
			% is a frame defined with its Z-axis as the axis of rotation of
			% this constraint.
			E_ac = eye(4);
			z = [0 0 1]';
			axis = this.bodyA.joint.axis;
			angle = acos(max(-1.0, min(axis'*z, 1.0)));
			E_ac(1:3,1:3) = se3.aaToMat(cross(axis,z),angle);
			E_ac(1:3,4) = this.xA;
			ntheta = 4;
			for k = 1 : ntheta
				theta = 2*pi*(k-1)/ntheta;
				s = sin(theta);
				c = cos(theta);
				nor_c = [c s 0]';
				pos_c = this.radius*nor_c;
				tan_c = [-s c 0]';
				for i = 1 : 2
					z = this.height*(i-1.5); % two ends of the cylinder
					pos_c(3) = z;
					bin_c = cross(tan_c,nor_c);
					this.contacts{end+1}.pos_a = E_ac(1:3,:)*[pos_c;1];
					this.contacts{end}.nor_a = -E_ac(1:3,1:3)*nor_c;
					this.contacts{end}.tan_a = E_ac(1:3,1:3)*tan_c;
					if i == 1
						bin_c = -bin_c;
					end
					this.contacts{end+1}.pos_a = E_ac(1:3,:)*[pos_c;1];
					this.contacts{end}.nor_a = E_ac(1:3,1:3)*bin_c;
					this.contacts{end}.tan_a = E_ac(1:3,1:3)*tan_c;
				end
			end
			% Index for tangent matrix
			nt = redmax.Scene.countT();
			t = length(this.contacts);
			this.idxT = nt + (1:t);
			nt = nt + t;
			redmax.Scene.countT(nt);
		end
		
		%%
		function computeContactMultiplier_(this,h,SPreg)
			% Computes the contact Lagrange multiplier
			nc = length(this.contacts);
			Minv = diag(1./[this.bodyA.I_i;this.bodyB.I_i]);
			N = zeros(nc,12);
			for i = 1 : nc
				contact = this.contacts{i};
				nor_a = contact.nor_a;
				pos_a = contact.pos_a;
				N(i,1:6) = nor_a'*se3.Gamma(pos_a);
				E_ba = this.bodyB.E_iw*this.bodyA.E_wi;
				nor_b = E_ba(1:3,1:3)*nor_a;
				pos_b = E_ba(1:3,:)*[pos_a;1];
				N(i,7:12) = -nor_b'*se3.Gamma(pos_b);
			end
			% Since fa = -N'*a/h, we need to multiply by h to get a
			H = N*Minv*N' + SPreg*eye(nc);
			H = 0.5*(H + H');
			f = (N*Minv*this.fcon)*h;
			a = H\f;
			for i = 1 : nc
				this.contacts{i}.a = a(i);
			end
		end
		
		%%
		function T = computeTangentMatrix_(this,T)
			% Assumes a revolute constraint!
			for i = 1 : length(this.contacts)
				contact = this.contacts{i};
				tan_a = contact.tan_a;
				pos_a = contact.pos_a;
				T(this.idxT(i),this.bodyA.idxM) = tan_a'*se3.Gamma(pos_a);
				E_ba = this.bodyB.E_iw*this.bodyA.E_wi;
				tan_b = E_ba(1:3,1:3)*tan_a;
				pos_b = E_ba(1:3,:)*[pos_a;1];
				T(this.idxT(i),this.bodyB.idxM) = -tan_b'*se3.Gamma(pos_b);
			end
		end
		
		%%
		function [bl,bu,idx] = computeFrictionLimits_(this,mu,SPathresh,bl,bu,idx)
			% Assumes a revolute constraint!
			% Combine normal and binormal, since they share the same
			% tangent.
			for i0 = 1 : 2 : length(this.contacts)
				i1 = i0 + 1;
				a0 = abs(this.contacts{i0}.a);
				a1 = abs(this.contacts{i1}.a);
				a = a0 + a1;
				% Use internal mu
				bl(this.idxT(i0)) = -mu(1)*a;
				bu(this.idxT(i0)) =  mu(1)*a;
				if a > SPathresh
					idx(end+1) = this.idxT(i0); %#ok<AGROW>
				end
			end
		end
		
		%%
		function draw_(this)
			% Assumes a revolute constraint!
			n = 8;
			E_ia = eye(4);
			z = [0 0 1]';
			axis = this.bodyA.joint.axis;
			angle = acos(max(-1.0, min(axis'*z, 1.0)));
			E_ia(1:3,1:3) = se3.aaToMat(cross(axis,z),angle);
			E_ia(1:3,4) = this.xA;
			E = this.bodyA.E_wi*E_ia;
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
			X = reshape(XYZ(1,:),2,n+1);
			Y = reshape(XYZ(2,:),2,n+1);
			Z = reshape(XYZ(3,:),2,n+1);
			surf(X,Y,Z,'FaceColor',this.bodyA.color);
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
			X = reshape(XYZ(1,:),2,n+1);
			Y = reshape(XYZ(2,:),2,n+1);
			Z = reshape(XYZ(3,:),2,n+1);
			surf(X,Y,Z,'FaceColor',this.bodyB.color);
			X = R.*cos(T);
			Y = R.*sin(T);
			Z = h/2*ones(size(R));
			X = reshape(X,1,2*(n+1));
			Y = reshape(Y,1,2*(n+1));
			Z = reshape(Z,1,2*(n+1));
			XYZ = E*[X;Y;Z;ones(1,2*(n+1))];
			X = reshape(XYZ(1,:),2,n+1);
			Y = reshape(XYZ(2,:),2,n+1);
			Z = reshape(XYZ(3,:),2,n+1);
			surf(X,Y,Z,'FaceColor',this.bodyB.color);
		end
	end
end
