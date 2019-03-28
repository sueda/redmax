classdef ConstraintFloor < redmax.Constraint
	% ConstraintFloor Inequality contact with the floor.
	% Assumptions for simplicity:
	% - The body is a sphere.
	% - There is one contact point per body.
	
	%%
	properties
		body % Body to apply the constraint on
		E    % The "floor" frame, with Z-up.
	end
	
	%%
	methods
		%%
		function this = ConstraintFloor(body)
			this = this@redmax.Constraint(0,0,1,0);
			this.name = [body.name,'-FLOOR'];
			this.body = body;
			this.E = eye(4);
		end
		
		%%
		function setTransform(this,E)
			this.E = E;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [Cm,Cmdot,cm,cmdot,cmddot] = computeJacIneqM_(this,Cm,Cmdot,cm,cmdot,cmddot)
			rows = this.idxIM;
			cols = this.body.idxM;
			this.idxQ = cols;
			% Assuming a sphere
			r = this.body.radius;
			x_w = this.body.E_wi(:,4); % body center in world frame
			x_f = this.E\x_w; % transformed to "floor" frame
			z = x_f(3);
			if z < r
				% Contact point is x_f, with radius subtracted from Z
				x = x_f;
				x(3) = x(3) - r;
				% Transform to world then to body
				x = this.body.E_wi\(this.E*x);
				G = se3.Gamma(x(1:3));
				nor = this.E(1:3,3);
				R = this.body.E_wi(1:3,1:3);
				Cm(rows,cols) = -nor'*R*G;
				cm(rows) = r - z;
				this.activeM = true;
			else
				this.activeM = false;
			end
		end
		
		%%
		function generateContactsCollision_(this)
			if ~this.activeM
				return;
			end
			% Assuming a sphere
			r = this.body.radius;
			E_wi = this.body.E_wi;
			R_wi = E_wi(1:3,1:3);
			E_wf = this.E;
			x_w = E_wi(:,4); % body center in world frame
			x_f = E_wf\x_w; % transformed to "floor" frame
			% Contact point is x_f, with radius subtracted from Z
			x = x_f;
			x(3) = x(3) - r;
			% Transform to world then to body
			x = E_wi\(E_wf*x);
			% Normal
			nor_w = E_wf(1:3,3);
			nor_a = R_wi'*nor_w;
			% Tangents (use X and Y from the floor)
			tan_w = E_wf(1:3,1:2);
			tan_a = R_wi'*tan_w;
			this.contacts{1}.pos_a = x(1:3);
			this.contacts{1}.nor_a = nor_a;
			this.contacts{1}.tan_a = tan_a; % two tangents
			% Index for tangent matrix
			nt = redmax.Scene.countT();
			this.idxT = nt + (1:2);
			nt = nt + 2;
			redmax.Scene.countT(nt);
		end
		
		%%
		function scatterForceIneqM_(this,Cmt,lm) %#ok<*INUSL>
			if ~this.activeM
				return;
			end
			this.contacts{1}.a = lm(this.idxIM);
		end
		
		%%
		function computeContactMultiplier_(this,h,SPreg) %#ok<INUSD>
			% Since the Lagrange multiplier passed in through the function
			% scatterForceIneqM_() was divided by h, we multiply by h here.
			% This is a bit of a hack!
			if ~this.activeM
				return;
			end
			this.contacts{1}.a = this.contacts{1}.a*h;
		end
		
		%%
		function T = computeTangentMatrix_(this,T)
			if ~this.activeM
				return;
			end
			contact = this.contacts{1};
			tan_a = contact.tan_a; % 3x2 matrix
			pos_a = contact.pos_a;
			T(this.idxT,this.body.idxM) = tan_a'*se3.Gamma(pos_a);
		end
		
		%%
		function [bl,bu,idx] = computeFrictionLimits_(this,mu,SPathresh,bl,bu,idx)
			if ~this.activeM
				return;
			end
			a = this.contacts{1}.a;
			% Use external mu
			bl(this.idxT(1)) = -mu(2)*a;
			bu(this.idxT(1)) =  mu(2)*a;
			bl(this.idxT(2)) = -mu(2)*a;
			bu(this.idxT(2)) =  mu(2)*a;
			if a > SPathresh
				idx(end+(1:2)) = this.idxT;
			end
		end
		
		%%
		function draw_(this)
			se3.drawAxis(this.E);
			s = 10;
			v = this.E(1:3,:)*[-s -s 0 1; s -s 0 1; s s 0 1; -s s 0 1]';
			f = [1 2 3 4];
			patch('Faces',f,'Vertices',v','FaceColor',[0.5 0.5 0.5])
		end
	end
end
