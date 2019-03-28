classdef SpringPointDirection < redmax.Spring
	%SpringPointDirection A linear force at a point in a particular direction
	
	%%
	properties
		stiffness
		body % Body to apply the wrench to
		x_i  % Application point in local coords
		d_w  % Force direction in world coords
	end
	
	methods
		function this = SpringPointDirection(body,x_i)
			this = this@redmax.Spring();
			this.body = body;
			this.x_i = x_i;
			this.d_w = [0 0 0]';
			this.stiffness = 0;
		end
		
		%%
		function setStiffness(this,stiffness)
			this.stiffness = stiffness;
		end
		
		%%
		function setDirection(this,d_w)
			this.d_w = d_w;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [f,K,D] = computeForceStiffnessDamping_(this,f,K,D)
			[f_,K_] = this.computeFK();
			idxM = this.body.idxM;
			f(idxM) = f(idxM) + f_;
			K(idxM,idxM) = K(idxM,idxM) + K_;
			% For recursive dynamics
			this.body.wext_i = this.body.wext_i + f_;
			this.body.Kmdiag = this.body.Kmdiag + K_;
		end
		
		%%
		function V = computeEnergy_(this,V)
			% TODO: this assumes that the point is the local origin
			%x_w = this.body.E_wi*[this.x_i;1]; % point in world coords
			%V = V - this.stiffness*x_w(1:3)'*this.d_w;
		end
		
		%%
		function y = computeStiffnessProd_(this,x,y)
			[~,K] = this.computeFK();
			idxM = this.body.idxM;
			y(idxM) = y(idxM) + K*x(idxM);
		end
		
		%%
		function draw_(this)
			x_w = this.body.E_wi*[this.x_i;1];
			xs = [x_w(1:3),x_w(1:3)+this.d_w];
			plot3(xs(1,:),xs(2,:),xs(3,:),'r-');
		end
	end
	
	%%
	methods (Access = private)
		%%
		function [f,K] = computeFK(this)
			f_w = this.stiffness*this.d_w;
			R = this.body.E_wi(1:3,1:3);
			G = se3.Gamma(this.x_i);
			f_i = R'*f_w;
			f = G'*f_i;
			% df = [ [xl]*[R'*fw]  0 ]
			%      [   [R'*fw]     0 ]
			% Not symmetric!
			df = zeros(6);
			fibrac = se3.brac(f_i);
			df(1:3,1:3) = se3.brac(this.x_i)*fibrac;
			df(4:6,1:3) = fibrac;
			K = 0.5*(df + df'); % symmetrize
			K = 0*K; % DISABLE stiffness matrix
		end
	end
end

