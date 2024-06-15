classdef TaskBDF1PointPos < redmax.TaskBDF1
	%TaskBDF1PointPos Sets a point on a rigid body to be at a target point.
	% The parameters are the fixed torques on all of the joints.
	
	%%
	properties
		t       % time at which the measurement will be taken
		body    % body to move to the target
		xlocal  % the local point on the body
		xtarget % the world target position
		pscale  % scale factor for the torques
		wpos    % position weight
	end
	
	%%
	methods
		%%
		function this = TaskBDF1PointPos(scene)
			nparams = 0;
			for i = 1 : length(scene.joints)
				nparams = nparams + scene.joints{i}.ndof;
			end
			this = this@redmax.TaskBDF1(scene,nparams);
		end
		
		%%
		function setTime(this,t)
			this.t = t;
		end
		
		%%
		function setBody(this,body)
			this.body = body;
		end
		
		%%
		function setPoint(this,xlocal)
			this.xlocal = xlocal;
		end
		
		%%
		function setTarget(this,xtarget)
			this.xtarget = xtarget;
		end
		
		%%
		function setScale(this,pscale)
			this.pscale = pscale;
		end
		
		%%
		function setWeights(this,wreg,wpos)
			this.wreg = wreg;
			this.wpos = wpos;
		end
		
		%% Applies the parameters to the current step
		function applyStep(this)
			% The indexing here assumes that we are controlling all of the
			% joints in the scene
			for i = 1 : length(this.scene.joints)
				this.scene.joints{i}.tau = this.pscale*this.p(this.scene.joints{i}.idxR);
			end
		end
		
		%% Computes objective value and derivatives after a sim step
		function calcStep(this)
			nm = redmax.Scene.countM();
			nr = redmax.Scene.countR();
			
			wp = this.wpos;
			
			k = this.scene.k;
			dt = this.t - this.scene.t;
			if abs(dt) < 1e-6
				% Objective
				q = this.scene.history(k).q;
				qdot = this.scene.history(k).qdot;
				this.scene.joints{1}.setQ(q,qdot);
				this.scene.joints{1}.update();
				E = this.body.E_wi;
				xworld = E(1:3,:)*[this.xlocal;1];
				dx = xworld - this.xtarget;
				this.P = this.P + wp*0.5*(dx'*dx);
				
				% Derivative
				% dPdq = dP/dx * dx/dqm * dqm/dqr
				%      = wp*dx * R*G * J
				dxdqm = zeros(3,nm);
				R = E(1:3,1:3);
				dxdqm(:,this.body.idxM) = R*se3.Gamma(this.xlocal);
				J = this.scene.history(k).J;

                nr = redmax.Scene.countR();
	            for i = 1 : nr
		            dJdqi = this.scene.history(k).dJdq(:,:,i);
                    J(:,i) =  J(:,i) + dJdqi*q;
                end

				this.dPdq{k} = J'*dxdqm'*dx*wp;
			else
				this.dPdq{k} = zeros(nr,1);
			end
			
			% dg/dp: the derivative of Newton function g wrt the parameters
			% With BDF1:
			%   g = M*(q1 - q0 - h*qdot0) - h^2*(fr + J'*fm + fqvv)
			% The parameters are the torques scaled by pscale, which go
			% into fr, so:
			%   dg/dp = -h^2*pscale*I
			h = this.scene.h;
			kk = (k-1)*nr + (1:nr); % indices for time step k
			this.dgdp(kk,:) = -h^2*this.pscale*speye(nr);
		end
		
		%% Draws the current task state
		function draw(this)
			x1 = this.xtarget;
			plot3(x1(1),x1(2),x1(3),'gs');
			xw = this.body.E_wi*[this.xlocal;1];
			xw1 = [xw(1:3),x1];
			plot3(xw1(1,:),xw1(2,:),xw1(3,:),'g');
		end
	end
	
	%%
	methods (Access = protected)
	end
end
