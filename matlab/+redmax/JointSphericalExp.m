classdef JointSphericalExp < redmax.Joint
	%JointSphericalExp A spherical joint parameterized with exponential
	%coordinates.
	%
	% Unlike Euler angles, exponential coordinates does not require
	% switching the coordinate chart. If we get close to a singularity, all
	% we need to do is to remap q to another value that produces the same
	% rotation matrix.
	
	%%
	properties
		radius        % Radius for contact
	end
	
	%%
	methods
		%%
		function this = JointSphericalExp(parent,body)
			this = this@redmax.Joint(parent,body,3);
			this.radius = 1.0;
		end
		
		%%
		function setGeometry(this,radius)
			this.radius = radius;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function reparam_(this)
			q0 = this.q(1:3);
			[q1,flag] = se3.reparam(q0);
			if flag
				%fprintf('%s reparam\n',this.name);
				qdot0 = this.qdot(1:3);
				S0 = this.S(1:3,1:3);
				% Update to compute new S (ignore new Sdot)
				this.q(1:3) = q1;
				this.update_();
				S1 = this.S(1:3,1:3);
				% Update again to compute new Sdot
				this.qdot(1:3) = S1\(S0*qdot0);
				this.update_();
			end
		end
		
		%%
		function update_(this)
			r = this.q(1:3);
			rdot = this.qdot(1:3);

			this.Q = eye(4);
			R = se3.exp(r);
			this.Q(1:3,1:3) = R;
			
			this.S = zeros(6,3);
			this.Sdot = zeros(6,3);
			rr = r'*r;
			if rr < 1e-9
				this.S(1:3,1:3) = eye(3);
			else
				% May want to remove some of the duplicate computations later
				Rdot = zeros(3);
				rbrac = se3.brac(r);
				IR = eye(3) - R;
				d = 1/rr;
				for i = 1 : 3
					Bi = r(i)*rbrac;
					Ci = se3.brac(rbrac*IR(:,i));
					Ai = (Bi + Ci)*d;
					dRdri = Ai*R;
					this.S(1:3,i) = se3.unbrac(R'*dRdri);
					Rdot = Rdot + dRdri*rdot(i);
				end
				rdotbrac = se3.brac(rdot);
				ddot = -2/(rr*rr)*(r'*rdot);
				for i = 1 : 3
					Bi = r(i)*rbrac;
					Ci = se3.brac(rbrac*IR(:,i));
					Ai = (Bi + Ci)*d;
					Bidot = rdot(i)*rbrac + r(i)*rdotbrac;
					Cidot = se3.brac(rdotbrac*IR(:,i) - rbrac*Rdot(:,i));
					Aidot = (Bidot + Cidot)*d + (Bi + Ci)*ddot;
					RdotAiR = Rdot'*Ai*R;
					this.Sdot(1:3,i) = se3.unbrac(RdotAiR + R'*Aidot*R - RdotAiR');
				end
			end
		end
		
		%%
		function generateContacts_(this)
			E_ij = this.body.E_ij;
			for k = 1 : 3
				R_jc = eye(3);
				if k == 2
					R_jc = R_jc(:,[2,3,1]);
				elseif k == 3
					R_jc = R_jc(:,[3,1,2]);
				end
				for i = 1 : 2
					R_jc = -R_jc; % negative and positive
					x = -this.radius*R_jc(:,1);
					E_jc = [R_jc,x; 0 0 0 1];
					E_ic = E_ij*E_jc;
					this.contacts{end+1}.pos_i = E_ic(1:3,4);
					this.contacts{end}.nor_i = E_ic(1:3,1);
					this.contacts{end}.tan_i = E_ic(1:3,2:3);
				end
			end
			% Index for tangent matrix
			nt = redmax.Scene.countT();
			t = 2*length(this.contacts);
			this.idxT = nt + (1:t);
			nt = nt + t;
			redmax.Scene.countT(nt);
		end
		
		%%
		function T = computeTangentMatrix_(this,T)
			for i = 1 : length(this.contacts)
				ii = 2*(i-1) + (1:2);
				contact = this.contacts{i};
				tan_i = contact.tan_i;
				pos_i = contact.pos_i;
				G_i = se3.Gamma(pos_i);
				T(this.idxT(ii(1)),this.body.idxM) = tan_i(:,1)'*G_i;
				T(this.idxT(ii(2)),this.body.idxM) = tan_i(:,2)'*G_i;
				if ~isempty(this.parent)
					E_pi = this.parent.body.E_iw*this.body.E_wi;
					tan_p = E_pi(1:3,:)*[tan_i;0 0];
					pos_p = E_pi(1:3,:)*[pos_i;1];
					G_p = se3.Gamma(pos_p);
					T(this.idxT(ii(1)),this.parent.body.idxM) = -tan_p(:,1)'*G_p;
					T(this.idxT(ii(2)),this.parent.body.idxM) = -tan_p(:,2)'*G_p;
				end
			end
		end
		
		%%
		function [bl,bu,idx] = computeFrictionLimits_(this,mu,SPathresh,bl,bu,idx)
			for i = 1 : length(this.contacts)
				ii = 2*(i-1) + (1:2);
				a = abs(this.contacts{i}.a);
				bl(this.idxT(ii(1))) = -mu*a;
				bu(this.idxT(ii(1))) =  mu*a;
				bl(this.idxT(ii(2))) = -mu*a;
				bu(this.idxT(ii(2))) =  mu*a;
				if a > SPathresh
					idx(end+(1:2)) = this.idxT(ii);
				end
			end
		end
		
		%%
		function draw_(this)
			E = this.E_wj;
			r = this.radius;
			S = diag([r,r,r,1]);
			n = 8;
			[X,Y,Z] = sphere(n);
			X = reshape(X,1,(n+1)^2);
			Y = reshape(Y,1,(n+1)^2);
			Z = reshape(Z,1,(n+1)^2);
			XYZ = E*S*[X;Y;Z;ones(1,(n+1)^2)];
			X = reshape(XYZ(1,:),n+1,n+1);
			Y = reshape(XYZ(2,:),n+1,n+1);
			Z = reshape(XYZ(3,:),n+1,n+1);
			surf(X,Y,Z,'FaceColor',this.body.color);
		end
	end
	
	%%
	methods (Static)
		%%
		function test()
			rng(0);
			body = redmax.BodyCuboid(1,[1 1 1]);
			joint = redmax.JointSphericalExp([],body);
			q0 = se3.reparam(randn(3,1));
			qdot0 = randn(3,1);
			joint.q(1:3) = q0;
			joint.qdot(1:3) = qdot0;
			joint.update_();
			Q = joint.Q;
			S = joint.S;
			Sdot = joint.Sdot;
			
			% Finite difference for S
			h = sqrt(eps);
			S_ = zeros(6,3);
			for i = 1 : 3
				joint.q = q0;
				joint.q(i) = joint.q(i) + h;
				joint.update_();
				S_(:,i) = se3.unbrac(joint.Q\(joint.Q - Q)/h);
			end
			norm(S - S_)
			
			% Finite difference for Sdot
			joint.q = q0 + h*qdot0;
			joint.update_();
			Sdot_ = (joint.S - S)/h;
			norm(Sdot - Sdot_)
			
			% How do we update the velocity if q is reparameterized?
			q0 = 2*pi*rand(3,1);
			joint.q(1:3) = q0;
			joint.qdot(1:3) = qdot0;
			joint.update_();
			Q0 = joint.Q;
			S0 = joint.S;
			Sdot0 = joint.Sdot;
			
			q1 = se3.reparam(q0);
			joint.q(1:3) = q1;
			joint.update_();
			Q1 = joint.Q;
			S1 = joint.S;
			Sdot1 = joint.Sdot;
			
			S0 = S0(1:3,1:3);
			S1 = S1(1:3,1:3);
			qdot1 = S1\(S0*qdot0);
		end
	end
end
