classdef SpringDamper < redmax.Spring
	%SpringDamper A spring-damper force between two points
	% The scalar force is
	%    f = k*(l-L)/L - d*v
	% where l is the current length, and v is current velocity.
	
	%%
	properties
		k     % Stiffness parameter
		d     % Damping parameter
		L     % Rest length
		body1 % Body to apply the wrench to (null implies world)
		body2 % Body to apply the wrench to (null implies world)
		x_1   % Application point in local coords
		x_2   % Application point in local coords
		xinit_1
		xinit_2
	end
	
	%%
	methods
		%%
		function this = SpringDamper(body1,x_1,body2,x_2)
			this = this@redmax.Spring();
			this.body1 = body1;
			this.body2 = body2;
			this.x_1 = x_1;
			this.x_2 = x_2;
			this.k = 1.0;
			this.d = 1.0;
			this.L = 0.0; % Use initial length if not set
			this.xinit_1 = x_1;
			this.xinit_2 = x_2;
		end
		
		%%
		function setStiffness(this,k)
			this.k = k;
		end
		
		%%
		function setDamping(this,d)
			this.d = d;
		end
		
		%%
		function setRestLength(this,L)
			this.L = L;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [f,K,D] = computeForceStiffnessDamping_(this,f,K,D)
			[f_,K_,D_] = this.computeFKD();
			idxM1 = [];
			idxM2 = [];
			if ~isempty(this.body1)
				idxM1 = this.body1.idxM;
				f1 = f_(1:6);
				K11 = K_(1:6,1:6);
				D11 = D_(1:6,1:6);
				f(idxM1) = f(idxM1) + f1;
				K(idxM1,idxM1) = K(idxM1,idxM1) + K11;
				D(idxM1,idxM1) = D(idxM1,idxM1) + D11;
				% For recursive dynamics
				this.body1.wext_i = this.body1.wext_i + f1;
				this.body1.Kmdiag = this.body1.Kmdiag + K11;
				this.body1.Dmdiag = this.body1.Dmdiag + D11;
			end
			if ~isempty(this.body2)
				idxM2 = this.body2.idxM;
				f2 = f_(7:12);
				K22 = K_(7:12,7:12);
				D22 = D_(7:12,7:12);
				f(idxM2) = f(idxM2) + f2;
				K(idxM2,idxM2) = K(idxM2,idxM2) + K22;
				D(idxM2,idxM2) = D(idxM2,idxM2) + D22;
				% For recursive dynamics
				this.body2.wext_i = this.body2.wext_i + f2;
				this.body2.Kmdiag = this.body2.Kmdiag + K22;
				this.body2.Dmdiag = this.body2.Dmdiag + D22;
			end
			if ~isempty(this.body1) && ~isempty(this.body2)
				K12 = K_(1:6,7:12);
				K21 = K_(7:12,1:6);
				D12 = D_(1:6,7:12);
				D21 = D_(7:12,1:6);
				K(idxM1,idxM2) = K(idxM1,idxM2) + K12;
				K(idxM2,idxM1) = K(idxM2,idxM1) + K21;
				D(idxM1,idxM2) = D(idxM1,idxM2) + D12;
				D(idxM2,idxM1) = D(idxM2,idxM1) + D21;
			end
		end
		
		%%
		function V = computeEnergy_(this,V)
			% Just the position part
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
			l = norm(x2_w - x1_w);
			if this.L == 0
				this.L = l;
			end
			V = V + 0.5*this.k*((l-this.L)/this.L)^2;
		end
		
		%%
		function y = computeStiffnessProd_(this,x,y)
			[~,K] = this.computeFKD();
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
		function y = computeDampingProd_(this,x,y)
			[~,~,D] = this.computeFKD();
			if ~isempty(this.body1)
				idxM1 = this.body1.idxM;
				y(idxM1) = y(idxM1) + D(1:6,1:6)*x(idxM1);
			end
			if ~isempty(this.body2)
				idxM2 = this.body2.idxM;
				y(idxM2) = y(idxM2) + D(7:12,7:12)*x(idxM2);
			end
			if ~isempty(this.body1) && ~isempty(this.body2)
				y(idxM1) = y(idxM1) + D(1:6,7:12)*x(idxM2);
				y(idxM2) = y(idxM2) + D(7:12,1:6)*x(idxM1);
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
			x1 = E1(1:3,:)*[x1;1];
			x2 = E2(1:3,:)*[x2;1];
			l = norm(x2 - x1);
			if this.L == 0
				this.L = l;
			end
			s = (l-this.L)/this.L;
			snorm = 0.5*(-s+1);
			cm = redmax.SpringMuscle.CM(1:end,:);
			cmax = size(cm,1);
			color = cm(max(1,min(round(snorm*cmax),cmax)),:);
			x12 = [x1 x2];
			plot3(x12(1,:),x12(2,:),x12(3,:),'-','Color',color,'LineWidth',3);
		end
	end
	
	%%
	methods (Access = private)
		%%
		function [f,K,D] = computeFKD(this)
			E1 = eye(4);
			E2 = eye(4);
			phi1 = zeros(6,1);
			phi2 = zeros(6,1);
			x1 = this.x_1;
			x2 = this.x_2;
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
			G1 = se3.Gamma(x1);
			G2 = se3.Gamma(x2);
			x1_w = E1(1:3,:)*[x1;1];
			x2_w = E2(1:3,:)*[x2;1];
			dx_w = x2_w - x1_w;
			l = norm(dx_w);
			if this.L == 0
				this.L = l;
			end
			v1_w = R1*G1*phi1;
			v2_w = R2*G2*phi2;
			v = (dx_w'/l)*(v2_w - v1_w);
			fx_1 = -G1'*R1'*dx_w;
			fx_2 =  G2'*R2'*dx_w;
			fx = [fx_1;fx_2];
			fn = (1/l)*fx;
			fs = this.k*(l-this.L)/this.L - this.d*v/this.L;
			f = -fs*fn;
			
			% Stiffness scalar part
			dfsdx1 = -this.k/this.L*dx_w'/l;
			dfsdE = [dfsdx1 * R1*G1, -dfsdx1 * R2*G2];
			K = fn*dfsdE;
			K = -0.5*(K + K'); % symmetrize
			
			% Damping scalar part
			dir_w = dx_w/l;
			dfmdv1 = this.d/this.L*dir_w';
			dfmdphi = [dfmdv1*R1*G1, -dfmdv1*R2*G2];
			D = -fn*dfmdphi;
		end
	end

end
