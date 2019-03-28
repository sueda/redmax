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
			this.axis = axis;
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
			this.Q = [se3.aaToMat(this.axis,this.q),[0 0 0]'; 0 0 0 1];
			this.S(1:3) = this.axis;
		end
		
		%%
		function generateContacts_(this)
			% The contact points will be generated in the 'A' frame, which
			% is a frame defined with its Z-axis as the axis of rotation of
			% this joint.
			E_ja = eye(4);
			z = [0 0 1]';
			angle = acos(max(-1.0, min(this.axis'*z, 1.0)));
			E_ja(1:3,1:3) = se3.aaToMat(cross(this.axis,z),angle);
			E_ia = this.body.E_ij * E_ja;
			E_ac = eye(4); % where is the contact frame wrt axis frame
			ntheta = 4;
			for k = 1 : ntheta
				theta = 2*pi*(k-1)/ntheta;
				s = sin(theta);
				c = cos(theta);
				nor = [c s 0]';
				pos = this.radius*nor;
				tan = [-s c 0]';
				for i = 1 : 2
					z = this.height*(i-1.5); % two ends of the cylinder
					pos(3) = z;
					bin = cross(tan,nor);
					E_ac(1:3,:) = [tan nor bin pos];
					E_ic = E_ia*E_ac;
					t = E_ic(1:3,1);
					n = E_ic(1:3,2);
					b = E_ic(1:3,3);
					x = E_ic(1:3,4);
					this.contacts{end+1}.pos_i = x;
					this.contacts{end}.nor_i = -n;
					this.contacts{end}.tan_i = t;
					if i == 1
						b = -b;
					end
					this.contacts{end+1}.pos_i = x;
					this.contacts{end}.nor_i = b;
					this.contacts{end}.tan_i = t;
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
		function T = computeTangentMatrix_(this,T)
			for i = 1 : length(this.contacts)
				contact = this.contacts{i};
				tan_i = contact.tan_i;
				pos_i = contact.pos_i;
				T(this.idxT(i),this.body.idxM) = tan_i'*se3.Gamma(pos_i);
				if ~isempty(this.parent)
					E_pi = this.parent.body.E_iw*this.body.E_wi;
					tan_p = E_pi(1:3,:)*[tan_i;0];
					pos_p = E_pi(1:3,:)*[pos_i;1];
					T(this.idxT(i),this.parent.body.idxM) = -tan_p'*se3.Gamma(pos_p);
				end
			end
		end
		
		%%
		function [bl,bu,idx] = computeFrictionLimits_(this,mu,SPathresh,bl,bu,idx)
			% Combine normal and binormal, since they share the same
			% tangent.
			for i0 = 1 : 2 : length(this.contacts)
				i1 = i0 + 1;
				a0 = abs(this.contacts{i0}.a);
				a1 = abs(this.contacts{i1}.a);
				a = a0 + a1;
				bl(this.idxT(i0)) = -mu*a;
				bu(this.idxT(i0)) =  mu*a;
				if a > SPathresh
					idx(end+1) = this.idxT(i0); %#ok<AGROW>
				end
			end
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
