classdef (Abstract) ForceSpringMultiPointGeneric < redmax.Force
	%ForceSpringMultiPointGeneric A generic spring routed along a sequence of points
	
	%%
	properties
		bodies % Bodies to apply the wrench to (null implies world)
		xls  % Application point in local coords
	end
	
	methods
		%%
		function this = ForceSpringMultiPointGeneric()
			this = this@redmax.Force();
			this.bodies = {};
			this.xls = {};
		end
		
		%%
		function addBodyPoint(this,body,xl)
			this.bodies{end+1} = body;
			this.xls{end+1} = xl;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [fr,fm,Kr,Km,Dr,Dm] = computeValues_(this,fr,fm,Kr,Km,Dr,Dm)
			% Compute quantities for each point
			npts = length(this.bodies);
			E = cell(1,npts);
			xl = cell(1,npts);
			G = cell(1,npts);
			phi = cell(1,npts);
			R = cell(1,npts);
			p = cell(1,npts);
			xw = cell(1,npts);
			vw = cell(1,npts);
			for k = 1 : npts
				E{k} = eye(4);
				phi{k} = zeros(6,1);
				xl{k} = this.xls{k};
				G{k} = se3.Gamma(xl{k});
				if ~isempty(this.bodies{k})
					% If a body is null, assume it is the world body
					E{k} = this.bodies{k}.E_wi;
					phi{k} = this.bodies{k}.phi;
				end
				R{k} = E{k}(1:3,1:3);
				p{k} = E{k}(1:3,4);
				xw{k} = R{k}*xl{k} + p{k};
				vw{k} = R{k}*G{k}*phi{k};
			end
			
			% Compute force vector
			fn = zeros(6*npts,1);
			l = 0;
			ldot = 0;
			for k = 1 : npts - 1
				idxL1 = (k-1)*6 + (1:6);
				idxL2 = (k-0)*6 + (1:6);
				idxL = [idxL1,idxL2];
				xw1 = xw{k};
				xw2 = xw{k+1};
				vw1 = vw{k};
				vw2 = vw{k+1};
				R1 = R{k};
				R2 = R{k+1};
				G1 = G{k};
				G2 = G{k+1};
				dx = xw2 - xw1;				
				dv = vw2 - vw1;				
				dxlen = norm(dx);
				l = l + dxlen;
				ldot = ldot + (dx'*dv)/dxlen;
				fx1 =  G1'*R1'*dx;
				fx2 = -G2'*R2'*dx;
				fx = [fx1;fx2];
				fn(idxL) = fn(idxL) + (1/dxlen)*fx;
			end
			[~,fs,dfsdl,dfsdldot] = this.computeSpringForce(l,ldot);
			f = fs*fn;
			
			% Copy to global
			for k1 = 1 : npts
				if ~isempty(this.bodies{k1})
					idxG1 = this.bodies{k1}.idxM;
					idxL1 = (k1-1)*6 + (1:6);
					fm(idxG1) = fm(idxG1) + f(idxL1);
				end
			end
			if nargout == 2
				return
			end
			
			% Compute stiffness and damping matrices
			I = eye(3);
			Kn = zeros(6*npts); % normalized vector part
			dfsdq = zeros(1,6*npts);
			dfsdqdot = zeros(1,6*npts);
			exb = se3.brac([1 0 0]');
			eyb = se3.brac([0 1 0]');
			ezb = se3.brac([0 0 1]');
			for k = 1 : npts - 1
				idxL1 = (k-1)*6 + (1:6);
				idxL2 = (k-0)*6 + (1:6);
				idxL = [idxL1,idxL2];
				xw1 = xw{k};
				xw2 = xw{k+1};
				vw1 = vw{k};
				vw2 = vw{k+1};
				R1 = R{k};
				R2 = R{k+1};
				G1 = G{k};
				G2 = G{k+1};
				phi1 = phi{k};
				phi2 = phi{k+1};
				dx = xw2 - xw1;				
				dv = vw2 - vw1;				
				dxlen = norm(dx);
				dxnor = dx/dxlen;
				
				% K scalar term
				dldq = dxnor'*[-R1*G1, R2*G2];
				dldotdq = ((I - dxnor*dxnor')/dxlen*dv)'*[-R1*G1, R2*G2];
				dldotdq([1,7]) = dldotdq([1,7]) + dxnor'*[-R1*exb*G1*phi1, R2*exb*G2*phi2]; % x rotation of q1 and q2
				dldotdq([2,8]) = dldotdq([2,8]) + dxnor'*[-R1*eyb*G1*phi1, R2*eyb*G2*phi2]; % y rotation of q1 and q2
				dldotdq([3,9]) = dldotdq([3,9]) + dxnor'*[-R1*ezb*G1*phi1, R2*ezb*G2*phi2]; % z rotation of q1 and q2
				dfsdq(idxL) = dfsdq(idxL) + dfsdl*dldq + dfsdldot*dldotdq;
				
				% K normalized vector term: K1
				fx1 =  G1'*R1'*dx;
				fx2 = -G2'*R2'*dx;
				fx = [fx1;fx2];
				d = -dx/dxlen^3;
				K1 = fx*[d'*R1*G1,-d'*R2*G2];
				
				% K normalized vector term: K2
				K2 = zeros(12);
				p1 = p{k};
				p2 = p{k+1};
				x1b = se3.brac(xl{k});
				x2b = se3.brac(xl{k+1});
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
				K2 = (1/dxlen)*K2;
				
				Kn(idxL,idxL) = Kn(idxL,idxL) + K1 + K2;
				
				% D
				d = dfsdldot*dxnor;
				dfsdqdot(idxL1) = dfsdqdot(idxL1) - d'*R1*G1;
				dfsdqdot(idxL2) = dfsdqdot(idxL2) + d'*R2*G2;
			end
			Ks = fn*dfsdq;
			K = Ks - fs*Kn;
			D = fn*dfsdqdot;
			
			% Copy to global
			for k1 = 1 : npts
				if ~isempty(this.bodies{k1})
					idxG1 = this.bodies{k1}.idxM;
					idxL1 = (k1-1)*6 + (1:6);
					Km(idxG1,idxG1) = Km(idxG1,idxG1) + K(idxL1,idxL1);
					Dm(idxG1,idxG1) = Dm(idxG1,idxG1) + D(idxL1,idxL1);
					for k2 = k1 + 1 : npts
						if ~isempty(this.bodies{k2})
							idxG2 = this.bodies{k2}.idxM;
							idxL2 = (k2-1)*6 + (1:6);
							Km(idxG1,idxG2) = Km(idxG1,idxG2) + K(idxL1,idxL2);
							Km(idxG2,idxG1) = Km(idxG2,idxG1) + K(idxL2,idxL1);
							Dm(idxG1,idxG2) = Dm(idxG1,idxG2) + D(idxL1,idxL2);
							Dm(idxG2,idxG1) = Dm(idxG2,idxG1) + D(idxL2,idxL1);
						end
					end
				end
			end
		end
		
		%%
		function V = computeEnergy_(this,V)
			% Compute quantities for each point
			npts = length(this.bodies);
			xw = cell(1,npts);
			vw = cell(1,npts);
			for k = 1 : npts
				E = eye(4);
				phi = zeros(6,1);
				xl = this.xls{k};
				G = se3.Gamma(xl);
				if ~isempty(this.bodies{k})
					% If a body is null, assume it is the world body
					E = this.bodies{k}.E_wi;
					phi = this.bodies{k}.phi;
				end
				R = E(1:3,1:3);
				p = E(1:3,4);
				xw{k} = R*xl + p;
				vw{k} = R*G*phi;
			end
			l = 0;
			ldot = 0;
			for k = 1 : npts - 1
				xw1 = xw{k};
				xw2 = xw{k+1};
				vw1 = vw{k};
				vw2 = vw{k+1};
				dx = xw2 - xw1;				
				dv = vw2 - vw1;				
				dxlen = norm(dx);
				l = l + dxlen;
				ldot = ldot + (dx'*dv)/dxlen;
			end
			V = V + this.computeSpringForce(l,ldot);
		end
		
		%%
		function draw_(this)
			npts = length(this.bodies);
			xw = zeros(3,npts);
			for k = 1 : npts
				E = eye(4);
				xl = this.xls{k};
				if ~isempty(this.bodies{k})
					% If a body is null, assume it is the world body
					E = this.bodies{k}.E_wi;
				end
				R = E(1:3,1:3);
				p = E(1:3,4);
				xw(:,k) = R*xl + p;
			end
			plot3(xw(1,:),xw(2,:),xw(3,:),'r-');
		end
	end
	
	methods (Abstract, Access = protected)
		%%
		[f,dfdl,dfdldot] = computeSpringForce(this,l,ldot);
	end
	
	methods (Static)
		%%
		function test()
			redmax.Scene.clear();
			density = 1;
			sides = [3 2 1];
			n = 3; % number of points (>= 2)
			cuboids = cell(1,n);
			force = redmax.ForceCable();
			force.setStiffness(1e3);
			force.setDamping(1e3);
			for k = 1 : n
				cuboids{k} = redmax.BodyCuboid(density,sides);
				cuboids{k}.idxM = (k-1)*6 + (1:6);
				xl = randn(3,1);
				phi = randn(6,1);
				cuboids{k}.E_wi = se3.randE();
				cuboids{k}.phi = phi;
				force.addBodyPoint(cuboids{k},xl);
			end
			force.computeRestLength();
			% Move the bodies around randomly
			for k = 1 : n
				cuboids{k}.E_wi = se3.randE();
			end
			redmax.Scene.countM(6*n);
			[~,fm,~,Km,~,Dm] = force.computeValues();
			% Finite difference test
			Km_ = zeros(6*n);
			for i = 1 : n
				E = cuboids{i}.E_wi;
				for j = 1 : 6
					ij = (i-1)*6 + j;
					dE = zeros(6,1);
					dE(j) = sqrt(eps);
					E_ = E*se3.exp(dE);
					cuboids{i}.E_wi = E_;
					[~,fm_] = force.computeValues();
					cuboids{i}.E_wi = E;
					Km_(:,ij) = (fm_ - fm)/sqrt(eps);
				end
			end
			redmax.Scene.printError('Km',Km_,Km);
			Dm_ = zeros(6*n);
			for i = 1 : n
				phi = cuboids{i}.phi;
				for j = 1 : 6
					ij = (i-1)*6 + j;
					phi_ = phi;
					phi_(j) = phi_(j) + sqrt(eps);
					cuboids{i}.phi = phi_;
					[~,fm_] = force.computeValues();
					cuboids{i}.phi = phi;
					Dm_(:,ij) = (fm_ - fm)/sqrt(eps);
				end
			end
			redmax.Scene.printError('Dm',Dm_,Dm);
		end
	end
end
