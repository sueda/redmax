classdef JointSplineCurve < redmax.Joint
	% Implementation of "Spline Joints for Multibody Dynamics" by Lee and
	% Terzopoulos 2008.
	% Uses cubic B-splines.
	% Assumes a cyclic spline curve.
	
	%%
	properties
		Cs
		dCs
	end
	
	%%
	properties (Constant)
		B = 1/6*[ % Bspline coeffs
			 1 -3  3 -1
			 4  0 -6  3
			 1  3  3 -3
			 0  0  0  1
			];
		B1 = 1/6*[ % Bspline coeffs, 1st column removed
			-3  3 -1
			 0 -6  3
			 3  3 -3
			 0  0  1
			];
		B2 = 1/6*[ % Bspline coeffs, 1st 2 columns removed
			 3 -1
			-6  3
			 3 -3
			 0  1
			];
	end
	
	%%
	methods
		%%
		function this = JointSplineCurve(parent,body)
			this = this@redmax.Joint(parent,body,1);
			this.Cs = {};
			this.dCs = zeros(6,0);
		end
		
		%%
		function addControlFrame(this,C)
			this.Cs{end+1} = C;
			ncfs = length(this.Cs);
			if ncfs >= 2
				C0 = this.Cs{ncfs-1};
				C1 = C;
				this.dCs(:,ncfs) = se3.unbrac(se3.log(C0\C1));
			end
			% Cyclic
			C0 = C;
			C1 = this.Cs{1};
			this.dCs(:,1) = se3.unbrac(se3.log(C0\C1));
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this)
			this.Q = this.evalQ(this.q(1));
			[this.S,dSdq] = this.evalS(this.q(1));
			this.Sdot = dSdq*this.qdot(1);
		end
		
		%%
		function draw_(this)
			ncfs = length(this.Cs);
			s = min(this.body.sides);
			if isempty(this.parent)
				E_wp = eye(4);
			else
				E_wp = this.parent.E_wj;
			end
			E_wj0 = E_wp*this.E_pj0;
			% Draw control frames
			nan3 = nan(3,1);
			xlines = zeros(3,0);
			ylines = zeros(3,0);
			zlines = zeros(3,0);
			for i = 1 : ncfs
				[xlines1,ylines1,zlines1] = se3.drawAxis(E_wj0*this.Cs{i},s);
				xlines = [xlines,xlines1,nan3]; %#ok<*AGROW>
				ylines = [ylines,ylines1,nan3];
				zlines = [zlines,zlines1,nan3];
			end
			plot3(xlines(1,:),xlines(2,:),xlines(3,:),'r','LineWidth',2);
			plot3(ylines(1,:),ylines(2,:),ylines(3,:),'g','LineWidth',2);
			plot3(zlines(1,:),zlines(2,:),zlines(3,:),'b','LineWidth',2);
			% Draw spline frames
			xlines = zeros(3,0);
			ylines = zeros(3,0);
			zlines = zeros(3,0);
			qmax = ncfs;
			for q = linspace(0,qmax,qmax*5)
				Q = evalQ(this,q);
				[xlines1,ylines1,zlines1] = se3.drawAxis(E_wj0*Q,s);
				xlines = [xlines,xlines1,nan3];
				ylines = [ylines,ylines1,nan3];
				zlines = [zlines,zlines1,nan3];
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
			ncfs = length(this.Cs);
			% Wrap around.
			qmax = ncfs;
			if q < 0
				q = q + qmax;
			elseif q >= qmax
				q = q - qmax;
			end
			k = floor(q); % starting control frame (0-index)
			if k >= ncfs
				k = k - 1; % overflow
			end
			q_ = q - k; % local q in [0,1]
			Q = this.Cs{k+1}; % convert to 1-index
			for i = 2 : 4
				ki = k + i;
				if ki > ncfs
					ki = ki - ncfs; % Wrap for cyclic
				end
				Bsum = redmax.JointSplineCurve.Bsum(i,q_);
				dC = se3.brac(this.dCs(:,ki));
				Q = Q*se3.exp(dC*Bsum);
			end
		end
		
		%%
		function [S,dSdq] = evalS(this,q)
			% Evaluates spline frame derivatives
			ncfs = length(this.Cs);
			% Wrap around.
			qmax = ncfs;
			if q < 0
				q = q + qmax;
			elseif q >= qmax
				q = q - qmax;
			end
			k = floor(q); % starting control frame (0-index)
			if k >= ncfs
				k = k - 1; % overflow
			end
			q_ = q - k; % local q in [0,1]
			for i = 2 : 4
				ki = k + i;
				if ki > ncfs
					ki = ki - ncfs; % Wrap for cyclic
				end
				dC = this.dCs(:,ki);
				dBsum = redmax.JointSplineCurve.dBsum(i,q_);
				d2Bsum = redmax.JointSplineCurve.d2Bsum(i,q_);
				if i == 2
					S = dC*dBsum;
					dSdq = dC*d2Bsum;
				else
					Bsum = redmax.JointSplineCurve.Bsum(i,q_);
					Ad = se3.Ad(se3.inv(se3.exp(dC*Bsum)));
					ad = se3.ad(S);
					S = dC*dBsum + Ad*S;
					dSdq = dC*d2Bsum + Ad*(dSdq + ad*dC*dBsum);
				end
			end
		end
	end

	%%
	methods (Static)
		%%
		function f = Bsum(i,q)
			% Evaluates Btilde
			B = redmax.JointSplineCurve.B;
			qvec = [1 q q^2 q^3]';
			if i == 1
				b = B(1,:) + B(2,:) + B(3,:) + B(4,:);
			elseif i == 2
				b = B(2,:) + B(3,:) + B(4,:);
			elseif i == 3
				b = B(3,:) + B(4,:);
			else
				b = B(4,:);
			end
			f = b*qvec;
		end
		
		%%
		function df = dBsum(i,q)
			% Evaluates dBtilde/dq
			B = redmax.JointSplineCurve.B1;
			qvec = [1 2*q 3*q^2]';
			if i == 1
				b = B(1,:) + B(2,:) + B(3,:) + B(4,:);
			elseif i == 2
				b = B(2,:) + B(3,:) + B(4,:);
			elseif i == 3
				b = B(3,:) + B(4,:);
			else
				b = B(4,:);
			end
			df = b*qvec;
		end
		
		%%
		function d2f = d2Bsum(i,q)
			% Evaluates d^2Btilde/dq^2
			B = redmax.JointSplineCurve.B2;
			qvec = [2 6*q]';
			if i == 1
				b = B(1,:) + B(2,:) + B(3,:) + B(4,:);
			elseif i == 2
				b = B(2,:) + B(3,:) + B(4,:);
			elseif i == 3
				b = B(3,:) + B(4,:);
			else
				b = B(4,:);
			end
			d2f = b*qvec;
		end

	end
end
