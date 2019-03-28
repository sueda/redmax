classdef JointSplineSurface < redmax.Joint
	% Implementation of "Spline Joints for Multibody Dynamics" by Lee and
	% Terzopoulos 2008.
	% Uses cubic B-spline tensor product patches.
	% Currently only works with a single patch of 4x4 control frames.
	
	%%
	properties
		cs
	end
	
	%%
	properties (Constant)
		B = 1/6*[ % Bspline coeffs
			 1 -3  3 -1
			 4  0 -6  3
			 1  3  3 -3
			 0  0  0  1
			];
		E = [ % 6 basis twists in Eq. (25)
			0     0     0     1     0     0
			0     0     0     0     1     0
			0     0     0     0     0     1
			1     0     0     0     0     0
			0     1     0     0     0     0
			0     0     1     0     0     0
			];
	end
	
	%%
	methods
		%%
		function this = JointSplineSurface(parent,body)
			this = this@redmax.Joint(parent,body,2);
			this.cs = zeros(4,4,6);
		end
		
		%%
		function setControlFrame(this,i,j,c)
			this.cs(i,j,:) = c;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this)
			this.Q = this.evalQ(this.q);
			[this.S,dSdq] = this.evalS(this.q);
			this.Sdot = dSdq(:,:,1)*this.qdot(1) + dSdq(:,:,2)*this.qdot(2);
		end
		
		%%
		function draw_(this)
			nan3 = nan(3,1);
			s = min(this.body.sides);
			if isempty(this.parent)
				E_wp = eye(4);
			else
				E_wp = this.parent.E_wj;
			end
			E_wj0 = E_wp*this.E_pj0;
			if false
				% Draw control frames
				xlines = zeros(3,0);
				ylines = zeros(3,0);
				zlines = zeros(3,0);
				for i1 = 1 : 4
					for i2 = 1 : 4
						q = this.cs(i1,i2,:);
						Q = eye(4);
						for k = 1 : 6
							Q = Q*se3.exp(this.E(:,k)*q(k));
						end
						[xlines1,ylines1,zlines1] = se3.drawAxis(E_wj0*Q,s);
						xlines = [xlines,xlines1,nan3]; %#ok<*AGROW>
						ylines = [ylines,ylines1,nan3];
						zlines = [zlines,zlines1,nan3];
					end
				end
				plot3(xlines(1,:),xlines(2,:),xlines(3,:),'r','LineWidth',2);
				plot3(ylines(1,:),ylines(2,:),ylines(3,:),'g','LineWidth',2);
				plot3(zlines(1,:),zlines(2,:),zlines(3,:),'b','LineWidth',2);
			end
			% Draw spline frames
			xlines = zeros(3,0);
			ylines = zeros(3,0);
			zlines = zeros(3,0);
			for q1 = linspace(0,1,5)
				for q2 = linspace(0,1,5)
					Q = evalQ(this,[q1,q2]');
					[xlines1,ylines1,zlines1] = se3.drawAxis(E_wj0*Q,s);
					xlines = [xlines,xlines1,nan3];
					ylines = [ylines,ylines1,nan3];
					zlines = [zlines,zlines1,nan3];
				end
			end
			% Draw current frame
			[xlines1,ylines1,zlines1] = se3.drawAxis(this.E_wj,2*s);
			xlines = [xlines,xlines1,nan3];
			ylines = [ylines,ylines1,nan3];
			zlines = [zlines,zlines1,nan3];
			plot3(xlines(1,:),xlines(2,:),xlines(3,:),'r','LineWidth',1);
			plot3(ylines(1,:),ylines(2,:),ylines(3,:),'g','LineWidth',1);
			plot3(zlines(1,:),zlines(2,:),zlines(3,:),'b','LineWidth',1);
		end
	end
	
	%%
	methods (Access = private)
		%%
		function Q = evalQ(this,q)
			% Evaluates spline frame
			Q = eye(4);
			for i = 1 : 6
				Ci = this.cs(:,:,i); % 4x4 control point values
				ei = redmax.JointSplineSurface.E(:,i);
				phi = redmax.JointSplineSurface.Cfun(Ci,q);
				Q = Q*se3.exp(ei*phi);
			end
		end
		
		%%
		function [S,dSdq] = evalS(this,q)
			% Evaluates spline frame derivatives
			S = zeros(6,2);
			dSdq = zeros(6,2,2);
			e1 = redmax.JointSplineSurface.E(:,1);
			C1 = this.cs(:,:,1);
			for i = 1 : 2
				dphi1i = redmax.JointSplineSurface.dCfun(C1,q,i);
				S(1:6,i) = e1*dphi1i;
				for j = 1 : 2
					d2phi1ij = redmax.JointSplineSurface.d2Cfun(C1,q,i,j);
					dSdq(1:6,i,j) = e1*d2phi1ij;
				end
			end
			for k = 2 : 6
				ek = redmax.JointSplineSurface.E(:,k);
				Ck = this.cs(:,:,k); % 4x4 control point values
				phik = redmax.JointSplineSurface.Cfun(Ck,q);
				Ad = se3.Ad(se3.inv(se3.exp(ek*phik)));
				for i = 1 : 2
					dphiki = redmax.JointSplineSurface.dCfun(Ck,q,i);
					ad = se3.ad(S(1:6,i));
					S(1:6,i) = ek*dphiki + Ad*S(1:6,i);
					for j = 1 : 2
						d2phikij = redmax.JointSplineSurface.d2Cfun(Ck,q,i,j);
						dphikj = redmax.JointSplineSurface.dCfun(Ck,q,j);
						dSdq(1:6,i,j) = ek*d2phikij + Ad*(dSdq(1:6,i,j) + ad*ek*dphikj);
					end
				end
			end
		end
	end

	%%
	methods (Static)
		%%
		function f = Cfun(C,q)
			B = redmax.JointSplineSurface.B;
			q1 = q(1);
			q2 = q(2);
			q1vec = [1 q1 q1^2 q1^3]';
			q2vec = [1 q2 q2^2 q2^3]';
			f = (q2vec'*B')*C*(B*q1vec);
		end
		
		%%
		function df = dCfun(C,q,i)
			B = redmax.JointSplineSurface.B;
			q1 = q(1);
			q2 = q(2);
			if i == 1
				q1vec = [0 1 2*q1 3*q1^2]';
				q2vec = [1 q2 q2^2 q2^3]';
			else
				q1vec = [1 q1 q1^2 q1^3]';
				q2vec = [0 1 2*q2 3*q2^2]';
			end
			df = (q2vec'*B')*C*(B*q1vec);
		end
		
		%%
		function d2f = d2Cfun(C,q,i,j)
			B = redmax.JointSplineSurface.B;
			q1 = q(1);
			q2 = q(2);
			if i == 1 && j == 1
				q1vec = [0 0 2 6*q1]';
				q2vec = [1 q2 q2^2 q2^3]';
			elseif i == 1 && j == 2
				q1vec = [0 1 2*q1 3*q1^2]';
				q2vec = [0 1 2*q2 3*q2^2]';
			elseif i == 2 && j == 1
				q1vec = [0 1 2*q1 3*q1^2]';
				q2vec = [0 1 2*q2 3*q2^2]';
			elseif i == 2 && j == 2
				q1vec = [1 q1 q1^2 q1^3]';
				q2vec = [0 0 2 6*q2]';
			end
			d2f = (q2vec'*B')*C*(B*q1vec);
		end

	end
end
