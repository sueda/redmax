classdef DeformableSpring < redmax.Deformable
	%DeformableSpring Non-zero rest-length serial spring
	%   A simple spring made up of a series of nodes. In the future, we may
	%   want to make this an abstract class and derive from it.
	
	%%
	properties
		nodes
		K
		body0
		body1
		r0
		r1
		mass
	end
	
	%%
	methods
		%%
		function this = DeformableSpring(nnodes)
			this = this@redmax.Deformable();
			global countD countCM;
			this.name = ['spring',num2str(countD)];
			cm = colormap('lines');
			this.color = cm(countCM,:);
			for i = 1 : nnodes
				this.nodes(i).x = zeros(3,1);
				this.nodes(i).v = zeros(3,1);
				this.nodes(i).a = zeros(3,1);
			end
			this.K = 0;
			this.mass = 1;
			countD = countD + 1;
			countCM = countCM + 1;
		end
		
		%%
		function setStiffness(this,K)
			this.K = K;
		end
		
		%%
		function setMass(this,mass)
			this.mass = mass;
		end
		
		%%
		function setAttachments(this,body0,r0,body1,r1)
			% Attaches this spring to body0 and body1.
			this.body0 = body0;
			this.body1 = body1;
			this.r0 = r0;
			this.r1 = r1;
		end
	end
	
	%%
	methods (Access = protected)
		
		%%
		function initGeometry_(this)
			% Sets the world positions of the nodes using the attachment
			% points. r0 and r1 are in local coords.
			if isempty(this.body0)
				E0 = eye(4);
			else
				E0 = this.body0.E_wi;
			end
			if isempty(this.body1)
				E1 = eye(4);
			else
				E1 = this.body1.E_wi;
			end
			x0 = E0*[this.r0;1];
			x1 = E1*[this.r1;1];
			x0 = x0(1:3);
			x1 = x1(1:3);
			% Set the nodal positions
			nnodes = length(this.nodes);
			for i = 1 : nnodes
				s = (i-1)/(nnodes-1);
				this.nodes(i).x = (1-s)*x0 + s*x1;
			end
			% Compute the rest lengths
			for i = 1 : nnodes-1
				x0 = this.nodes(i).x;
				x1 = this.nodes(i+1).x;
				dx = x1 - x0;
				this.nodes(i).L = norm(dx); % rest length
			end
		end
		
		%%
		function countDofs_(this)
			% Counts maximal and reduced DOFs
			% For non-rigid DOFs, we need both maximal and reduced DOFs,
			% and the Jacobian must pass them through with the identity
			% matrices.
			nm = redmax.Scene.countM();
			nr = redmax.Scene.countR();
			for i = 1 : length(this.nodes)
				this.nodes(i).idxM = nm + (1:3);
				this.nodes(i).idxR = nr + (1:3);
				nm = nm + 3;
				nr = nr + 3;
			end
			redmax.Scene.countM(nm);
			redmax.Scene.countR(nr);
		end
		
		%%
		function y = gatherDofs_(this,y)
			% Gathers q and qdot into y
			nr = redmax.Scene.countR();
			for i = 1 : length(this.nodes)
				idxR = this.nodes(i).idxR;
				y(idxR) = this.nodes(i).x;
				y(nr + idxR) = this.nodes(i).v;
			end
		end
		
		%%
		function ydot = gatherDDofs_(this,ydot)
			% Gathers qdot and qddot into ydot
			nr = redmax.Scene.countR();
			for i = 1 : length(this.nodes)
				idxR = this.nodes(i).idxR;
				ydot(idxR) = this.nodes(i).v;
				ydot(nr + idxR) = this.nodes(i).a;
			end
		end
		
		%%
		function scatterDofs_(this,y)
			% Scatters q and qdot from y
			nr = redmax.Scene.countR();
			for i = 1 : length(this.nodes)
				idxR = this.nodes(i).idxR;
				this.nodes(i).x(1:3) = y(idxR);
				this.nodes(i).v(1:3) = y(nr + idxR);
			end
		end
		
		%%
		function scatterDDofs_(this,ydot)
			% Scatters qdot and qddot from ydot
			nr = redmax.Scene.countR();
			for i = 1 : length(this.nodes)
				idxR = this.nodes(i).idxR;
				this.nodes(i).v(1:3) = ydot(idxR);
				this.nodes(i).a(1:3) = ydot(nr + idxR);
			end
		end
		
		%%
		function [J,Jdot] = computeJacobian_(this,J,Jdot)
			for i = 1 : length(this.nodes)
				J(this.nodes(i).idxM,this.nodes(i).idxR) = eye(3);
			end
		end
		
		%%
		function [M,f] = computeMassGrav_(this,grav,M,f)
			% Computes maximal mass matrix and force vector
			nnodes = length(this.nodes);
			m = this.mass/nnodes;
			I = eye(3);
			for i = 1 : nnodes
				rows = this.nodes(i).idxM;
				M(rows,rows) = m*I;
				f(rows) = f(rows) + m*grav;
			end
			for i = 1 : nnodes-1
				rows0 = this.nodes(i).idxM;
				rows1 = this.nodes(i+1).idxM;
				x0 = this.nodes(i).x;
				x1 = this.nodes(i+1).x;
				dx = x1 - x0;
				l = sqrt(dx'*dx);
				L = this.nodes(i).L;
				e = (l - L)/L;
				fs = this.K*e*(1/L)*dx/l;
				f(rows0) = f(rows0) + fs;
				f(rows1) = f(rows1) - fs;
			end
		end
		
		%%
		function [f,D] = computeForceDamping_(this,f,D)
			% Computes maximal damping vector and matrix
			for i = 1 : length(this.nodes)
				rows = this.nodes(i).idxM;
				f(rows) = f(rows) - this.damping*this.nodes(i).v;
				D(rows,rows) = D(rows,rows) + this.damping*eye(3);
			end
		end
		
		%%
		function [T,V] = computeEnergies_(this,grav,T,V)
			nnodes = length(this.nodes);
			m = this.mass/nnodes;
			for i = 1 : nnodes
				x = this.nodes(i).x;
				v = this.nodes(i).v;
				T = T + 0.5*m*(v'*v);
				V = V - m*grav'*x;
			end
			for i = 1 : nnodes-1
				x0 = this.nodes(i).x;
				x1 = this.nodes(i+1).x;
				L = this.nodes(i).L;
				dx = x1 - x0;
				l = sqrt(dx'*dx);
				e = (l - L)/L;
				V = V + 0.5*this.K*e^2;
			end
		end
		
		%%
		function draw_(this)
			xs = [this.nodes.x];
			plot3(xs(1,:),xs(2,:),xs(3,:),'o-','Color',this.color);
		end
	end
end

