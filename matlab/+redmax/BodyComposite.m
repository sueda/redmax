classdef BodyComposite < redmax.Body
	%BodyComposite A body made of multiple bodies
	
	%%
	properties
		bodies
		Es     % Transforms of the bodies
		I_c    % Composite inertia
	end
	
	%%
	methods
		%%
		function this = BodyComposite(density)
			this = this@redmax.Body(density);
			this.bodies = {};
			this.Es = {};
		end
		
		%%
		function addBody(this,body,E)
			% E is the transform of body wrt parent joint.
			this.bodies{end+1} = body;
			this.Es{end+1} = E;
		end
		
		%%
		function E_jc = computeInertiaFrame(this)
			% Computes the composite inertia and frame.
			% The computations will be performed wrt the parent joint's
			% frame. We can't use world coords because we don't know where
			% these bodies are going to go.
			n = length(this.bodies);
			% Total mass
			mtotal = 0;
			for k = 1 : n
				body = this.bodies{k};
				body.computeInertia_();
				mtotal = mtotal + body.I_i(6);
			end
			% Center of mass
			xc = zeros(3,1);
			for k = 1 : n
				body = this.bodies{k};
				E = this.Es{k};
				xc = xc + (body.I_i(6)/mtotal)*E(1:3,4);
			end
			% Construct a frame at xc
			E_jc = eye(4);
			E_jc(1:3,4) = xc;
			% Sum the inertias at E_wc
			this.I_c = zeros(6);
			for k = 1 : n
				body = this.bodies{k};
				E_jk = this.Es{k};
				E_kj = se3.inv(E_jk);
				E_kc = E_kj * E_jc;
				Ad_kc = se3.Ad(E_kc);
				this.I_c = this.I_c + Ad_kc'*diag(body.I_i)*Ad_kc;
			end
			% Compute the rotated frame
			J = this.I_c(1:3,1:3);
			[V,D] = eig(J);
			if cross(V(:,1),V(:,2))'*V(:,3) < 0
				% Switch to right-handed coord frame
				V(1:3,3) = -V(1:3,3);
			end
			E_jc(1:3,1:3) = V;
			this.I_c = [diag(D)' mtotal mtotal mtotal]';
			% Transform Es to be the transform of each body wrt the
			% composite inertia
			E_cj = se3.inv(E_jc);
			for k = 1 : n
				E_jb = this.Es{k};
				E_cb = E_cj * E_jb;
				this.Es{k} = E_cb;
			end
		end
		
		%%
		function computeInertia_(this)
			% Computes inertia at body
			this.I_i = this.I_c;
		end
		
		%%
		function draw_(this)
			se3.drawAxis(this.E_wi,5);
			for k = 1 : length(this.bodies)
				body = this.bodies{k};
				body.E_wi = this.E_wi * this.Es{k};
				body.draw_();
			end
		end
		
		%%
		function s = getAxisSize(this)
			s = 0;
			for k = 1 : length(this.bodies)
				body = this.bodies{k};
				s = max(s,body.getAxisSize());
			end
		end
	end
end
