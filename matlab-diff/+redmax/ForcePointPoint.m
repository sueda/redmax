classdef ForcePointPoint < redmax.Force
	%ForcePointPoint A linear zero rest-length force between two points
	
	%%
	properties
		stiffness
		damping
		body1 % Body to apply the wrench to (null implies world)
		body2 % Body to apply the wrench to (null implies world)
		x_1  % Application point in local coords
		x_2  % Application point in local coords
	end
	
	methods
		function this = ForcePointPoint(body1,x_1,body2,x_2)
			this = this@redmax.Force();
			if ~isempty(body1)
				this.name = body1.name;
			else
				this.name = 'NULL';
			end
			this.name = [this.name,'-'];
			if ~isempty(body2)
				this.name = [this.name,body2.name];
			else
				this.name = [this.name,'NULL'];
			end
			this.body1 = body1;
			this.body2 = body2;
			this.x_1 = x_1;
			this.x_2 = x_2;
			this.stiffness = 1;
			this.damping = 0;
		end
		
		%%
		function setStiffness(this,stiffness)
			this.stiffness = stiffness;
		end
		
		%%
		function setDamping(this,damping)
			this.damping = damping;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [fr,fm,Kr,Km,Dr,Dm] = computeValues_(this,fr,fm,Kr,Km,Dr,Dm)
			if ~isempty(this.body1)
				E1 = this.body1.E_wi;
				R1 = E1(1:3,1:3);
				p1 = E1(1:3,4);
				xl1 = this.x_1; % local position to apply the force to
				xw1 = R1*xl1 + p1; % point transformed to world
				G1 = se3.Gamma(xl1);
				vl1 = G1*this.body1.phi; % local velocity
				vw1 = R1*vl1; % velocity transformed to world
			else
				xw1 = this.x_1;
				vw1 = zeros(3,1);
			end
			if ~isempty(this.body2)
				E2 = this.body2.E_wi;
				R2 = E2(1:3,1:3);
				p2 = E2(1:3,4);
				xl2 = this.x_2; % local position to apply the force to
				xw2 = R2*xl2 + p2; % point transformed to world
				G2 = se3.Gamma(xl2);
				vl2 = G2*this.body2.phi; % local velocity
				vw2 = R2*vl2; % velocity transformed to world
			else
				xw2 = this.x_2;
				vw2 = zeros(3,1);
			end
			dx = xw2 - xw1;
			dv = vw2 - vw1;
			I = eye(3);
			Z = zeros(3);
			ks = this.stiffness;
			kd = this.damping;
			f = (ks*dx + kd*dv);
			% Force vector
			if ~isempty(this.body1)
				idx1 = this.body1.idxM;
				fm(idx1) = fm(idx1) + G1'*R1'*f;
			end
			if ~isempty(this.body2)
				idx2 = this.body2.idxM;
				fm(idx2) = fm(idx2) - G2'*R2'*f;
			end
			if nargout == 2
				return
			end
			% Stiffness and damping matrices
			if ~isempty(this.body1)
				Km(idx1,idx1) = Km(idx1,idx1) + ks*G1'*[se3.brac(R1'*(xw2 - p1)), -I];
				Km(idx1,idx1) = Km(idx1,idx1) + kd*G1'*[se3.brac(R1'*vw2), Z];
				Dm(idx1,idx1) = Dm(idx1,idx1) - kd*(G1'*G1);
			end
			if ~isempty(this.body2)
				Km(idx2,idx2) = Km(idx2,idx2) + ks*G2'*[se3.brac(R2'*(xw1 - p2)), -I];
				Km(idx2,idx2) = Km(idx2,idx2) + kd*G2'*[se3.brac(R2'*vw1), Z];
				Dm(idx2,idx2) = Dm(idx2,idx2) - kd*(G2'*G2);
			end
			if ~isempty(this.body1) && ~isempty(this.body2)
				Km(idx1,idx2) = Km(idx1,idx2) + ks*G1'*R1'*R2*[-se3.brac(xl2), I];
				Km(idx2,idx1) = Km(idx2,idx1) + ks*G2'*R2'*R1*[-se3.brac(xl1), I];
				Km(idx1,idx2) = Km(idx1,idx2) - kd*G1'*R1'*R2*[se3.brac(vl2), Z];
				Km(idx2,idx1) = Km(idx2,idx1) - kd*G2'*R2'*R1*[se3.brac(vl1), Z];
				Dm(idx1,idx2) = Dm(idx1,idx2) + kd*G1'*R1'*R2*G2;
				Dm(idx2,idx1) = Dm(idx2,idx1) + kd*G2'*R2'*R1*G1;
			end
		end
		
		%%
		function V = computeEnergy_(this,V)
			E1 = eye(4);
			E2 = eye(4);
			x1 = this.x_1;
			x2 = this.x_2;
			if ~isempty(this.body1)
				E1 = this.body1.E_wi;
			end
			if ~isempty(this.body2)
				E2 = this.body2.E_wi;
			end
			x1_w = E1(1:3,:)*[x1;1];
			x2_w = E2(1:3,:)*[x2;1];
			dx_w = x2_w - x1_w;
			V = V + 0.5*this.stiffness*(dx_w'*dx_w);
		end
		
		%%
		function draw_(this)
			E1 = eye(4);
			E2 = eye(4);
			x1 = this.x_1;
			x2 = this.x_2;
			if ~isempty(this.body1)
				E1 = this.body1.E_wi;
			end
			if ~isempty(this.body2)
				E2 = this.body2.E_wi;
			end
			x1_w = E1(1:3,:)*[x1;1];
			x2_w = E2(1:3,:)*[x2;1];
			xs = [x1_w x2_w];
			plot3(xs(1,:),xs(2,:),xs(3,:),'r-');
		end
	end
	
	methods (Static)
		%%
		function test()
			redmax.Scene.clear();
			density = 1;
			sides = [3 2 1];
			cuboid1 = redmax.BodyCuboid(density,sides);
			cuboid2 = redmax.BodyCuboid(density,sides);
			cuboid1.idxM = 1:6;
			cuboid2.idxM = 7:12;
			redmax.Scene.countM(12);
			E1 = se3.randE();
			E2 = se3.randE();
			phi1 = randn(6,1);
			phi2 = randn(6,1);
			x1 = randn(3,1);
			x2 = randn(3,1);
			cuboid1.E_wi = E1;
			cuboid2.E_wi = E2;
			cuboid1.phi = phi1;
			cuboid2.phi = phi2;
			force = redmax.ForcePointPoint(cuboid1,x1,cuboid2,x2);
			force.setStiffness(0e3);
			force.setDamping(1e2);
			[~,fm,~,Km,~,Dm] = force.computeValues();
			% Finite difference test
			Km_ = zeros(12);
			Dm_ = zeros(12);
			for i = 1 : 12
				% K
				phi = zeros(6,1);
				if i <= 6
					phi(i) = 1;
					E1_ = E1*se3.exp(sqrt(eps)*phi);
					cuboid1.E_wi = E1_;
				else
					phi(i-6) = 1;
					E2_ = E2*se3.exp(sqrt(eps)*phi);
					cuboid2.E_wi = E2_;
				end
				[~,fm_] = force.computeValues();
				cuboid1.E_wi = E1;
				cuboid2.E_wi = E2;
				Km_(:,i) = (fm_ - fm)/sqrt(eps);
				% D
				if i <= 6
					phi1_ = phi1;
					phi1_(i) = phi1_(i) + sqrt(eps);
					cuboid1.phi = phi1_;
				else
					phi2_ = phi2;
					phi2_(i-6) = phi2_(i-6) + sqrt(eps);
					cuboid2.phi = phi2_;
				end
				[~,fm_] = force.computeValues();
				cuboid1.phi = phi1;
				cuboid2.phi = phi2;
				Dm_(:,i) = (fm_ - fm)/sqrt(eps);
			end
			redmax.Scene.printError('Km',Km_,Km);
			redmax.Scene.printError('Dm',Dm_,Dm);
		end
	end
end
