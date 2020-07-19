classdef ForceGroundCuboid < redmax.Force
	% ForceGroundCuboid Contact with the ground.
	% Currently only works with cuboids.
	
	%%
	properties
		cuboid % BodyCuboid to apply the constraint on
		E    % The "ground" frame, with Z-up.
		kn   % normal stiffness
		kt   % tangential stiffness
		mu   % coefficient of friction
	end
	
	%%
	methods
		%%
		function this = ForceGroundCuboid(cuboid)
			this = this@redmax.Force();
			this.name = [cuboid.name,'-GROUND'];
			this.cuboid = cuboid;
			this.E = eye(4);
			this.kn = 1;
			this.kt = 0;
			this.mu = 0;
		end
		
		%%
		function setTransform(this,E)
			this.E = E;
		end
		
		%%
		function setStiffness(this,kn,kt)
			this.kn = kn;
			this.kt = kt;
		end
		
		%%
		function setFriction(this,mu)
			this.mu = mu;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [fr,fm,Kr,Km,Dr,Dm] = computeValues_(this,fr,fm,Kr,Km,Dr,Dm)
			idxM = this.cuboid.idxM;
			xg = this.E(1:3,4); % ground origin
			ng = this.E(1:3,3); % ground normal
			N = ng*ng';
			I = eye(3);
			Z = zeros(3);
			T = I - N;
			R = this.cuboid.E_wi(1:3,1:3); % body rotation
			p = this.cuboid.E_wi(1:3,4); % body position
			phi = this.cuboid.phi;
			e1b = se3.brac([1 0 0]');
			e2b = se3.brac([0 1 0]');
			e3b = se3.brac([0 0 1]');
			RNR = R'*N*R;
			pxgtmp = se3.brac(R'*N*(p - xg));
			% Collision detection
			S = eye(4);
			S(1:3,1:3) = diag(0.5*this.cuboid.sides);
			xl = S*[
				-1 -1 -1  1
				-1 -1  1  1
				-1  1 -1  1
				-1  1  1  1
				 1 -1 -1  1
				 1 -1  1  1
				 1  1 -1  1
				 1  1  1  1
				]';
			xw = this.cuboid.E_wi*xl;
			for i = 1 : 8
				xli = xl(1:3,i);
				xwi = xw(1:3,i);
				RNRxl = RNR*xli;
				xlbrac = se3.brac(xli);
				% penetration depth
				d = ng'*(xwi - xg);
				if d > 0
					% No collision
					continue
				end
				
				% Contact force
				fc = -this.kn*ng*d;
				G = se3.Gamma(xli);
				fm(idxM) = fm(idxM) + G'*R'*fc;
				if nargout > 2
					% Contact stiffness
					tmp = -[e1b*RNRxl,e2b*RNRxl,e3b*RNRxl] - RNR*xlbrac + pxgtmp;
					Km(idxM,idxM) = Km(idxM,idxM) - this.kn*G'*[tmp RNR];
				end
				
				% Friction force
				if this.mu == 0
					continue;
				end
				xwdot = R*G*phi;
				a = T*xwdot;
				anorm = norm(a);
				if this.mu*abs(this.kn*d) > this.kt*anorm
					% Static friction force
					fs = -this.kt*a;
					fm(idxM) = fm(idxM) + G'*R'*fs;
					if nargout > 2
						% Static friction damping and stiffness
						D = -this.kt*G'*R'*T*R*G;
						Gphi = G*phi;
						B = R'*T*R;
						K = -this.kt*G'*([(B*e1b-e1b*B)*Gphi, (B*e2b-e2b*B)*Gphi, (B*e3b-e3b*B)*Gphi, Z]);
						Dm(idxM,idxM) = Dm(idxM,idxM) + D;
						Km(idxM,idxM) = Km(idxM,idxM) + K;
					end
				else
					% dynamic friction force
					mukn = this.mu*this.kn;
					t = a/anorm;
					fd = -mukn*d*t;
					fm(idxM) = fm(idxM) + G'*R'*fd;
					if nargout > 2
						% dynamic friction damping and stiffness
						A = (a'*a*I - a*a')/norm(a)^3;
						D = -mukn*G'*R'*d*A*T*R*G;
						K1 = -d*[e1b*R'*t, e2b*R'*t, e3b*R'*t, Z];
						K2 = R'*t*ng'*R*G;
						K3 = -d*R'*A*T*R*[se3.brac(G*phi), Z];
						K = -mukn*G'*(K1 + K2 + K3);
						Dm(idxM,idxM) = Dm(idxM,idxM) + D;
						Km(idxM,idxM) = Km(idxM,idxM) + K;
					end
				end
			end
		end
		
		%%
		function V = computeEnergy_(this,V)
			xg = this.E(1:3,4); % ground origin
			ng = this.E(1:3,3); % ground normal
			% Collision detection
			S = eye(4);
			S(1:3,1:3) = diag(0.5*this.cuboid.sides);
			xl = S*[
				-1 -1 -1  1
				-1 -1  1  1
				-1  1 -1  1
				-1  1  1  1
				 1 -1 -1  1
				 1 -1  1  1
				 1  1 -1  1
				 1  1  1  1
				]';
			xw = this.cuboid.E_wi*xl;
			for i = 1 : 8
				x = xw(1:3,i);
				% penetration depth
				d = ng'*(x - xg);
				if d > 0
					% No collision
					continue
				end
				V = V + 0.5*this.kn*(d'*d);
			end
		end
		
		%%
		function draw_(this)
			se3.drawAxis(this.E);
			s = 10;
			v = this.E(1:3,:)*[-s -s 0 1; s -s 0 1; s s 0 1; -s s 0 1]';
			f = [1 2 3 4];
			patch('Faces',f,'Vertices',v','FaceColor',[0.5 0.5 0.5])
		end
	end
	
	methods (Static)
		%%
		function test()
			% Contact force
			redmax.Scene.clear();
			density = 1;
			sides = [3 2 1];
			cuboid = redmax.BodyCuboid(density,sides);
			cuboid.idxM = 1:6;
			redmax.Scene.countM(6);
			force = redmax.ForceGroundCuboid(cuboid);
			kn = 1e3;
			h = 1e-2;
			force.setStiffness(kn,h*kn);
			force.setFriction(0.5);
			E = se3.randE();
			cuboid.E_wi = E;
			cuboid.phi = randn(6,1);
			[~,fm,~,Km] = force.computeValues();
			% Finite difference test
			Km_ = zeros(6);
			for i = 1 : 6
				phi = zeros(6,1);
				phi(i) = 1;
				E_ = E*se3.exp(sqrt(eps)*phi);
				cuboid.E_wi = E_;
				[~,fm_] = force.computeValues();
				Km_(:,i) = (fm_ - fm)/sqrt(eps);
			end
			redmax.Scene.printError('Km',Km_,Km);
		end
		
		%%
		function test2()
			% Test derivative of R'*N*R
			E = se3.randE();
			R = E(1:3,1:3);
			n = randn(3,1);
			n = n/norm(n);
			x = randn(3,1);
			N = n*n';
			RNR = R'*N*R;
			f = RNR*x;
			e1b = se3.brac([1 0 0]');
			e2b = se3.brac([0 1 0]');
			e3b = se3.brac([0 0 1]');
			K = -[e1b*f, e2b*f, e3b*f] - RNR*se3.brac(x);
			% Finite difference test
			K_ = zeros(3,3);
			for i = 1 : 3
				w = zeros(3,1);
				w(i) = 1;
				R_ = R*expm(se3.brac(sqrt(eps)*w));
				f_ = R_'*N*R_*x;
				K_(:,i) = (f_ - f)/sqrt(eps);
			end
			redmax.Scene.printError('K',K_,K);
		end
		
		%%
		function test3()
			% Static friction
			E = se3.randE();
			E(1:3,1:3) = eye(3);
			phi = randn(6,1);
			kt = 1e3;
			R = E(1:3,1:3);
			ng = randn(3,1);
			ng = ng/norm(ng);
			xl = randn(3,1);
			N = ng*ng';
			T = eye(3) - N;
			G = se3.Gamma(xl);
			f = -kt*G'*R'*T*R*G*phi;
			D = -kt*G'*R'*T*R*G;
			Gphi = G*phi;
			B = R'*T*R;
			e1b = se3.brac([1 0 0]');
			e2b = se3.brac([0 1 0]');
			e3b = se3.brac([0 0 1]');
			Z = zeros(3,3);
			K = -kt*G'*([(B*e1b-e1b*B)*Gphi, (B*e2b-e2b*B)*Gphi, (B*e3b-e3b*B)*Gphi, Z]);
			% Finite difference test
			D_ = zeros(6,6);
			K_ = zeros(6,6);
			for i = 1 : 6
				phi_ = phi;
				phi_(i) = phi_(i) + sqrt(eps);
				f_ = -kt*G'*R'*T*R*G*phi_;
				D_(:,i) = (f_ - f)/sqrt(eps);
				phi_ = zeros(6,1);
				phi_(i) = sqrt(eps);
				E_ = E*se3.exp(se3.brac(phi_));
				R_ = E_(1:3,1:3);
				f_ = -kt*G'*R_'*T*R_*G*phi;
				K_(:,i) = (f_ - f)/sqrt(eps);
			end
			redmax.Scene.printError('D',D_,D);
			redmax.Scene.printError('K',K_,K);
		end
		
		%%
		function test4()
			% dynamic friction
			mu = 0.8;
			E = se3.randE();
			E(1:3,1:3) = eye(3);
			phi = randn(6,1);
			kt = 1e3;
			R = E(1:3,1:3);
			p = E(1:3,4);
			ng = randn(3,1);
			ng = ng/norm(ng);
			xg = randn(3,1);
			xl = randn(3,1);
			xw = R*xl + p;
			N = ng*ng';
			I = eye(3);
			T = I - N;
			G = se3.Gamma(xl);
			xwdot = R*G*phi;
			a = T*xwdot;
			t = a/norm(a);
			d = ng'*(xw - xg);
			f = -mu*kt*G'*R'*d*t;
			A = (a'*a*I - a*a')/norm(a)^3;
			D = -mu*kt*G'*R'*d*A*T*R*G;
			e1b = se3.brac([1 0 0]');
			e2b = se3.brac([0 1 0]');
			e3b = se3.brac([0 0 1]');
			Z = zeros(3,3);
			K1 = -d*[e1b*R'*t, e2b*R'*t, e3b*R'*t, Z];
			K2 = R'*t*ng'*R*G;
			K3 = -d*R'*A*T*R*[se3.brac(G*phi), Z];
			K = -mu*kt*G'*(K1 + K2 + K3);
			% Finite difference test
			D_ = zeros(6,6);
			K_ = zeros(6,6);
			for i = 1 : 6
				phi_ = phi;
				phi_(i) = phi_(i) + sqrt(eps);
				xwdot_ = R*G*phi_;
				a_ = T*xwdot_;
				t_ = a_/norm(a_);
				f_ = -mu*kt*G'*R'*d*t_;
				D_(:,i) = (f_ - f)/sqrt(eps);
				phi_ = zeros(6,1);
				phi_(i) = sqrt(eps);
				E_ = E*se3.exp(se3.brac(phi_));
				R_ = E_(1:3,1:3);
				p_ = E_(1:3,4);
				xwdot_ = R_*G*phi;
				a_ = T*xwdot_;
				t_ = a_/norm(a_);
				f_ = -mu*kt*G'*R_'*(ng'*(R_*xl + p_ - xg))*t_;
				K_(:,i) = (f_ - f)/sqrt(eps);
			end
			redmax.Scene.printError('D',D_,D);
			redmax.Scene.printError('K',K_,K);
		end
		
		%%
		function test5
			% dynamic friction individual terms
			E = se3.randE();
			E(1:3,1:3) = eye(3);
			phi = randn(6,1);
			R = E(1:3,1:3);
			p = E(1:3,4);
			ng = randn(3,1);
			ng = ng/norm(ng);
			xg = randn(3,1);
			xl = randn(3,1);
			xw = R*xl + p;
			N = ng*ng';
			I = eye(3);
			T = I - N;
			G = se3.Gamma(xl);
			xwdot = R*G*phi;
			a = T*xwdot;
			A = (a'*a*I - a*a')/norm(a)^3;
			t = a/norm(a);
			d = ng'*(xw - xg);
			f = d*R'*t;
			e1b = se3.brac([1 0 0]');
			e2b = se3.brac([0 1 0]');
			e3b = se3.brac([0 0 1]');
			Z = zeros(3,3);
			K1 = -d*[e1b*R'*t, e2b*R'*t, e3b*R'*t, Z];
			K2 = R'*t*ng'*R*G;
			K3 = -d*R'*A*T*R*[se3.brac(G*phi), Z];
			% Finite difference test
			K1_ = zeros(3,6);
			K2_ = zeros(3,6);
			K3_ = zeros(3,6);
			for i = 1 : 6
				phi_ = zeros(6,1);
				phi_(i) = sqrt(eps);
				E_ = E*se3.exp(se3.brac(phi_));
				R_ = E_(1:3,1:3);
				p_ = E_(1:3,4);
				xw_ = R_*xl + p_;
				d_ = ng'*(xw_ - xg);
				xwdot_ = R_*G*phi;
				a_ = T*xwdot_;
				t_ = a_/norm(a_);
				f1_ = d*R_'*t;
				f2_ = d_*R'*t;
				f3_ = d*R'*t_;
				K1_(:,i) = (f1_ - f)/sqrt(eps);
				K2_(:,i) = (f2_ - f)/sqrt(eps);
				K3_(:,i) = (f3_ - f)/sqrt(eps);
			end
			redmax.Scene.printError('K1',K1_,K1);
			redmax.Scene.printError('K2',K2_,K2);
			redmax.Scene.printError('K3',K3_,K3);
		end
	end
end
