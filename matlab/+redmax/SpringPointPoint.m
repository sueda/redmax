classdef SpringPointPoint < redmax.Spring
	%SpringPointPoint A linear force between two points
	
	%%
	properties
		stiffness
		body1 % Body to apply the wrench to (null implies world)
		body2 % Body to apply the wrench to (null implies world)
		x_1  % Application point in local coords
		x_2  % Application point in local coords
	end
	
	methods
		function this = SpringPointPoint(body1,x_1,body2,x_2)
			this = this@redmax.Spring();
			this.body1 = body1;
			this.body2 = body2;
			this.x_1 = x_1;
			this.x_2 = x_2;
			this.stiffness = 0;
		end
		
		%%
		function setStiffness(this,stiffness)
			this.stiffness = stiffness;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [f,K,D] = computeForceStiffnessDamping_(this,f,K,D)
			[f_,K_] = this.computeFK();
			idxM1 = [];
			idxM2 = [];
			if ~isempty(this.body1)
				idxM1 = this.body1.idxM;
				f1 = f_(1:6);
				K11 = K_(1:6,1:6);
				f(idxM1) = f(idxM1) + f1;
				K(idxM1,idxM1) = K(idxM1,idxM1) + K11;
				% For recursive dynamics
				this.body1.wext_i = this.body1.wext_i + f1;
				this.body1.Kmdiag = this.body1.Kmdiag + K11;
			end
			if ~isempty(this.body2)
				idxM2 = this.body2.idxM;
				f2 = f_(7:12);
				K22 = K_(7:12,7:12);
				f(idxM2) = f(idxM2) + f2;
				K(idxM2,idxM2) = K(idxM2,idxM2) + K22;
				% For recursive dynamics
				this.body2.wext_i = this.body2.wext_i + f2;
				this.body2.Kmdiag = this.body2.Kmdiag + K22;
			end
			if ~isempty(this.body1) && ~isempty(this.body2)
				K12 = K_(1:6,7:12);
				K21 = K_(7:12,1:6);
				K(idxM1,idxM2) = K(idxM1,idxM2) + K12;
				K(idxM2,idxM1) = K(idxM2,idxM1) + K21;
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
			f_w = x2_w - x1_w;
			V = V + 0.5*this.stiffness*(f_w'*f_w);
		end
		
		%%
		function y = computeStiffnessProd_(this,x,y)
			[~,K] = this.computeFK();
			if ~isempty(this.body1)
				idxM1 = this.body1.idxM;
				y(idxM1) = y(idxM1) + K(1:6,1:6)*x(idxM1);
			end
			if ~isempty(this.body2)
				idxM2 = this.body2.idxM;
				y(idxM2) = y(idxM2) + K(7:12,7:12)*x(idxM2);
			end
			if ~isempty(this.body1) && ~isempty(this.body2)
				y(idxM1) = y(idxM1) + K(1:6,7:12)*x(idxM2);
				y(idxM2) = y(idxM2) + K(7:12,1:6)*x(idxM1);
			end
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
	
	%%
	methods (Access = private)
		%%
		function [f,K] = computeFK(this)
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
			R1 = E1(1:3,1:3);
			R2 = E2(1:3,1:3);
			p1 = E1(1:3,4);
			p2 = E2(1:3,4);
			G1 = se3.Gamma(x1);
			G2 = se3.Gamma(x2);
			x1_w = E1(1:3,:)*[x1;1];
			x2_w = E2(1:3,:)*[x2;1];
			f_w = this.stiffness*(x2_w - x1_w);
			w_1 =  G1'*R1'*f_w;
			w_2 = -G2'*R2'*f_w;
			f = [w_1;w_2];
			% The Hessian matrix is not symmetric! (see checkDerivatives.m)
			dw = zeros(12);
			x1b = se3.brac(x1);
			x2b = se3.brac(x2);
			I = eye(3);
			R2R1 = R2'*R1;
			R1R2 = R2R1';
			dw( 4: 6, 1: 3) = se3.brac(R1'*(p1 - x2_w));
			dw( 1: 3, 1: 3) = x1b*dw( 4: 6, 1: 3);
			dw(10:12, 1: 3) = R2R1*x1b;
			dw( 7: 9, 1: 3) = x2b*dw(10:12, 1: 3);
			dw( 4: 6, 4: 6) = I;
			dw( 1: 3, 4: 6) = x1b;
			dw(10:12, 4: 6) = -R2R1;
			dw( 7: 9, 4: 6) = x2b*dw(10:12, 4: 6);
			dw( 4: 6, 7: 9) = R1R2*x2b;
			dw( 1: 3, 7: 9) = x1b*dw( 4: 6, 7: 9);
			dw(10:12, 7: 9) = se3.brac(R2'*(p2 - x1_w));
			dw( 7: 9, 7: 9) = x2b*dw(10:12, 7: 9);
			dw( 4: 6,10:12) = -R1R2;
			dw( 1: 3,10:12) = x1b*dw( 4: 6,10:12);
			dw(10:12,10:12) = I;
			dw( 7: 9,10:12) = x2b;
			K = -0.5*this.stiffness*(dw + dw'); % symmetrize
		end
	end
end

