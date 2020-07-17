classdef ForcePointPoint < redmax.Force
	%ForcePointPoint A linear force between two points
	
	%%
	properties
		stiffness
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
		end
		
		%%
		function setStiffness(this,stiffness)
			this.stiffness = stiffness;
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
				xl1 = this.x_1;
				xl1brac = se3.brac(xl1);
				xw1 = R1*xl1 + p1; % point transformed to world
			else
				xw1 = this.x_1;
			end
			if ~isempty(this.body2)
				E2 = this.body2.E_wi;
				R2 = E2(1:3,1:3);
				p2 = E2(1:3,4);
				xl2 = this.x_2; % local position to apply the force to
				xl2brac = se3.brac(xl2);
				xw2 = R2*xl2 + p2; % point transformed to world
			else
				xw2 = this.x_2;
			end
			dx = xw2 - xw1;
			I = eye(3);
			k = this.stiffness;
			if ~isempty(this.body1)
				idxM1 = this.body1.idxM;
				fm(idxM1(1:3)) = fm(idxM1(1:3)) + k*xl1brac*R1'*dx;
				fm(idxM1(4:6)) = fm(idxM1(4:6)) + k*R1'*dx;
			end
			if ~isempty(this.body2)
				idxM2 = this.body2.idxM;
				fm(idxM2(1:3)) = fm(idxM2(1:3)) - k*xl2brac*R2'*dx;
				fm(idxM2(4:6)) = fm(idxM2(4:6)) - k*R2'*dx;
			end
			if nargout == 2
				return
			end
			if ~isempty(this.body1)
				idxM1 = this.body1.idxM;
				% 1,1
				Km(idxM1(1:3),idxM1(1:3)) = Km(idxM1(1:3),idxM1(1:3)) - k*xl1brac*se3.brac(R1'*(p1 - xw2));
				Km(idxM1(1:3),idxM1(4:6)) = Km(idxM1(1:3),idxM1(4:6)) - k*xl1brac;
				Km(idxM1(4:6),idxM1(1:3)) = Km(idxM1(4:6),idxM1(1:3)) - k*se3.brac(R1'*(p1 - xw2));
				Km(idxM1(4:6),idxM1(4:6)) = Km(idxM1(4:6),idxM1(4:6)) - k*I;
			end
			if ~isempty(this.body2)
				idxM2 = this.body2.idxM;
				% 2,2
				Km(idxM2(1:3),idxM2(1:3)) = Km(idxM2(1:3),idxM2(1:3)) - k*xl2brac*se3.brac(R2'*(p2 - xw1));
				Km(idxM2(1:3),idxM2(4:6)) = Km(idxM2(1:3),idxM2(4:6)) - k*xl2brac;
				Km(idxM2(4:6),idxM2(1:3)) = Km(idxM2(4:6),idxM2(1:3)) - k*se3.brac(R2'*(p2 - xw1));
				Km(idxM2(4:6),idxM2(4:6)) = Km(idxM2(4:6),idxM2(4:6)) - k*I;
			end
			if ~isempty(this.body1) && ~isempty(this.body2)
				% 1,2
				Km(idxM1(1:3),idxM2(1:3)) = Km(idxM1(1:3),idxM2(1:3)) - k*xl1brac*R1'*R2*xl2brac;
				Km(idxM1(1:3),idxM2(4:6)) = Km(idxM1(1:3),idxM2(4:6)) + k*xl1brac*R1'*R2;
				Km(idxM1(4:6),idxM2(1:3)) = Km(idxM1(4:6),idxM2(1:3)) - k*R1'*R2*xl2brac;
				Km(idxM1(4:6),idxM2(4:6)) = Km(idxM1(4:6),idxM2(4:6)) + k*R1'*R2;
				% 2,1
				Km(idxM2(1:3),idxM1(1:3)) = Km(idxM2(1:3),idxM1(1:3)) - k*xl2brac*R2'*R1*xl1brac;
				Km(idxM2(1:3),idxM1(4:6)) = Km(idxM2(1:3),idxM1(4:6)) + k*xl2brac*R2'*R1;
				Km(idxM2(4:6),idxM1(1:3)) = Km(idxM2(4:6),idxM1(1:3)) - k*R2'*R1*xl1brac;
				Km(idxM2(4:6),idxM1(4:6)) = Km(idxM2(4:6),idxM1(4:6)) + k*R2'*R1;
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
			x1 = randn(3,1);
			x2 = randn(3,1);
			force = redmax.ForcePointPoint(cuboid1,x1,cuboid2,x2);
			force.setStiffness(1e3);
			E1 = se3.randE();
			E2 = se3.randE();
			cuboid1.E_wi = E1;
			cuboid2.E_wi = E2;
			[~,fm,~,Km] = force.computeValues();
			% Finite difference test
			Km_ = zeros(12);
			for i = 1 : 12
				phi = zeros(6,1);
				if i <= 6
					phi(i) = 1;
					E1_ = E1*se3.exp(sqrt(eps)*phi);
					cuboid1.E_wi = E1_;
					cuboid2.E_wi = E2;
				else
					phi(i-6) = 1;
					E2_ = E2*se3.exp(sqrt(eps)*phi);
					cuboid1.E_wi = E1;
					cuboid2.E_wi = E2_;
				end
				[~,fm_] = force.computeValues();
				Km_(:,i) = (fm_ - fm)/sqrt(eps);
			end
			redmax.Scene.printError('Km',Km_,Km);
		end
	end
end
