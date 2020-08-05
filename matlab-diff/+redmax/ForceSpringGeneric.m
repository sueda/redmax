classdef (Abstract) ForceSpringGeneric < redmax.Force
	%ForceSpringGeneric A generic spring between two points
	
	%%
	properties
		body1 % Body to apply the wrench to (null implies world)
		body2 % Body to apply the wrench to (null implies world)
		x_1  % Application point in local coords
		x_2  % Application point in local coords
	end
	
	methods
		%%
		function this = ForceSpringGeneric(body1,x_1,body2,x_2)
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
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [fr,fm,Kr,Km,Dr,Dm] = computeValues_(this,fr,fm,Kr,Km,Dr,Dm)
			E1 = eye(4);
			E2 = eye(4);
			xl1 = this.x_1;
			xl2 = this.x_2;
			G1 = se3.Gamma(this.x_1);
			G2 = se3.Gamma(this.x_2);
			phi1 = zeros(6,1);
			phi2 = zeros(6,1);
			if ~isempty(this.body1)
				E1 = this.body1.E_wi;
				phi1 = this.body1.phi;
				idxM1 = this.body1.idxM;
			end
			if ~isempty(this.body2)
				E2 = this.body2.E_wi;
				phi2 = this.body2.phi;
				idxM2 = this.body2.idxM;
			end
			R1 = E1(1:3,1:3);
			R2 = E2(1:3,1:3);
			p1 = E1(1:3,4);
			p2 = E2(1:3,4);
			xw1 = R1*xl1 + p1;
			xw2 = R2*xl2 + p2;
			vw1 = R1*G1*phi1;
			vw2 = R2*G2*phi2;
			dx = xw2 - xw1;
			l = norm(dx);
			dv = vw2 - vw1;
			ldot = (dx'*dv)/l;
			[~,fs,dfsdl,dfsdldot] = this.computeSpringForce(l,ldot);
			
			% force vector
			fx_1 =  G1'*R1'*dx;
			fx_2 = -G2'*R2'*dx;
			f = (fs/l)*[fx_1;fx_2]; % generalized force
			if ~isempty(this.body1)
				fm(idxM1) = fm(idxM1) + f(1:6);
			end
			if ~isempty(this.body2)
				fm(idxM2) = fm(idxM2) + f(7:12);
			end
			if nargout == 2
				return
			end
			
			% stiffness and damping matrices
			
			% K1
			I = eye(3);
			dldq = dx'/l*[-R1*G1, R2*G2];
			exb = se3.brac([1 0 0]');
			eyb = se3.brac([0 1 0]');
			ezb = se3.brac([0 0 1]');
			dldotdq = ((dx'*dx*I - dx*dx')/l^3*dv)'*[-R1*G1, R2*G2];
			dldotdq([1,7]) = dldotdq([1,7]) + dx'/l*[-R1*exb*G1*phi1, R2*exb*G2*phi2]; % x rotation of q1 and q2
			dldotdq([2,8]) = dldotdq([2,8]) + dx'/l*[-R1*eyb*G1*phi1, R2*eyb*G2*phi2]; % y rotation of q1 and q2
			dldotdq([3,9]) = dldotdq([3,9]) + dx'/l*[-R1*ezb*G1*phi1, R2*ezb*G2*phi2]; % z rotation of q1 and q2
			dfsdq = dfsdl*dldq + dfsdldot*dldotdq;
			K1 = [fx_1;fx_2]*(dfsdq/l - fs/l^2*dldq);
			
			% K2
			K2 = zeros(12);
			x1b = se3.brac(xl1);
			x2b = se3.brac(xl2);
			R2R1 = R2'*R1;
			R1R2 = R2R1';
			K2( 4: 6, 1: 3) = se3.brac(R1'*(p1 - xw2));
			K2( 1: 3, 1: 3) = x1b*K2( 4: 6, 1: 3);
			K2(10:12, 1: 3) = R2R1*x1b;
			K2( 7: 9, 1: 3) = x2b*K2(10:12, 1: 3);
			K2( 4: 6, 4: 6) = I;
			K2( 1: 3, 4: 6) = x1b;
			K2(10:12, 4: 6) = -R2R1;
			K2( 7: 9, 4: 6) = x2b*K2(10:12, 4: 6);
			K2( 4: 6, 7: 9) = R1R2*x2b;
			K2( 1: 3, 7: 9) = x1b*K2( 4: 6, 7: 9);
			K2(10:12, 7: 9) = se3.brac(R2'*(p2 - xw1));
			K2( 7: 9, 7: 9) = x2b*K2(10:12, 7: 9);
			K2( 4: 6,10:12) = -R1R2;
			K2( 1: 3,10:12) = x1b*K2( 4: 6,10:12);
			K2(10:12,10:12) = I;
			K2( 7: 9,10:12) = x2b;
			K2 = -(fs/l)*K2;
			
			K = K1 + K2;
			
			% D
			d_w = dfsdldot*dx/l^2;
			D = -[fx_1;fx_2]*[d_w'*R1*G1,-d_w'*R2*G2];

			if ~isempty(this.body1)
				Km(idxM1,idxM1) = Km(idxM1,idxM1) + K(1:6,1:6);
				Dm(idxM1,idxM1) = Dm(idxM1,idxM1) + D(1:6,1:6);
			end
			if ~isempty(this.body2)
				Km(idxM2,idxM2) = Km(idxM2,idxM2) + K(7:12,7:12);
				Dm(idxM2,idxM2) = Dm(idxM2,idxM2) + D(7:12,7:12);
			end
			if ~isempty(this.body1) && ~isempty(this.body2)
				Km(idxM1,idxM2) = Km(idxM1,idxM2) + K(1:6,7:12);
				Km(idxM2,idxM1) = Km(idxM2,idxM1) + K(7:12,1:6);
				Dm(idxM1,idxM2) = Dm(idxM1,idxM2) + D(1:6,7:12);
				Dm(idxM2,idxM1) = Dm(idxM2,idxM1) + D(7:12,1:6);
			end
		end
		
		%%
		function V = computeEnergy_(this,V)
			E1 = eye(4);
			E2 = eye(4);
			xl1 = this.x_1;
			xl2 = this.x_2;
			G1 = se3.Gamma(this.x_1);
			G2 = se3.Gamma(this.x_2);
			phi1 = zeros(6,1);
			phi2 = zeros(6,1);
			if ~isempty(this.body1)
				E1 = this.body1.E_wi;
				phi1 = this.body1.phi;
			end
			if ~isempty(this.body2)
				E2 = this.body2.E_wi;
				phi2 = this.body2.phi;
			end
			R1 = E1(1:3,1:3);
			R2 = E2(1:3,1:3);
			p1 = E1(1:3,4);
			p2 = E2(1:3,4);
			xw1 = R1*xl1 + p1;
			xw2 = R2*xl2 + p2;
			vw1 = R1*G1*phi1;
			vw2 = R2*G2*phi2;
			dx = xw2 - xw1;
			l = norm(dx);
			dv = vw2 - vw1;
			ldot = (dx'*dv)/l;
			V = V + this.computeSpringForce(l,ldot);
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
			xw1 = E1(1:3,:)*[x1;1];
			xw2 = E2(1:3,:)*[x2;1];
			xs = [xw1 xw2];
			plot3(xs(1,:),xs(2,:),xs(3,:),'r-');
		end
	end
	
	methods (Abstract, Access = protected)
		%%
		[V,f,dfdl,dfdldot] = computeSpringForce(this,l,ldot);
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
			force = redmax.ForceSpringDamped(cuboid1,x1,cuboid2,x2);
			force.setStiffness(1e3);
			force.setDamping(1e3);
			E1 = se3.randE();
			E2 = se3.randE();
			phi1 = randn(6,1);
			phi2 = randn(6,1);
			cuboid1.E_wi = E1;
			cuboid2.E_wi = E2;
			cuboid1.phi = phi1;
			cuboid2.phi = phi2;
			[~,fm,~,Km,~,Dm] = force.computeValues();
			% Finite difference test
			Km_ = zeros(12);
			for i = 1 : 12
				dE = zeros(6,1);
				if i <= 6
					dE(i) = sqrt(eps);
					E1_ = E1*se3.exp(dE);
					cuboid1.E_wi = E1_;
				else
					dE(i-6) = sqrt(eps);
					E2_ = E2*se3.exp(dE);
					cuboid2.E_wi = E2_;
				end
				[~,fm_] = force.computeValues();
				cuboid1.E_wi = E1;
				cuboid2.E_wi = E2;
				Km_(:,i) = (fm_ - fm)/sqrt(eps);
			end
			redmax.Scene.printError('Km',Km_,Km);
			Dm_ = zeros(12);
			for i = 1 : 12
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
			redmax.Scene.printError('Dm',Dm_,Dm);
		end
	end
end
