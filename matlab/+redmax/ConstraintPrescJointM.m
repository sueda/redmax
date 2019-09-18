classdef ConstraintPrescJointM < redmax.Constraint
	% ConstraintPrescJointM Prescribes body motion in maximal coords
	% Only works with revolute joints and with velocity-level Euler
	
	%%
	properties
		joint
		vel   % velocity-level integrator?
		q
		qdot
		qddot
	end
	
	%%
	methods
		%%
		function this = ConstraintPrescJointM(joint,vel)
			this = this@redmax.Constraint(1,0,0,0);
			this.name = [joint.name,'-presc'];
			this.joint = joint;
			this.vel = vel;
			this.q = 0;
			this.qdot = 0;
			this.qddot = 0;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [Gm,Gmdot,gm,gmdot,gmddot] = computeJacEqM_(this,Gm,Gmdot,gm,gmdot,gmddot)
			% Assumes both bodies exist
			rows = this.idxEM;
			colsI = this.joint.body.idxM;
			colsP = this.joint.parent.body.idxM;
			this.idxQ = [colsI colsP];
			Ad_ji = this.joint.body.Ad_ji;
			Ad_iw = this.joint.body.Ad_iw;
			Ad_wp = this.joint.parent.body.Ad_wi;
			GmI = Ad_ji;
			GmP = -Ad_ji*Ad_iw*Ad_wp;
			if sum(this.joint.axis == [1 0 0]') == 3
				% Revolute around X-axis
				Gm(rows,colsI) = GmI(1,:);
				Gm(rows,colsP) = GmP(1,:);
			elseif sum(this.joint.axis == [0 1 0]') == 3
				% Revolute around Y-axis
				Gm(rows,colsI) = GmI(2,:);
				Gm(rows,colsP) = GmP(2,:);
			elseif sum(this.joint.axis == [0 0 1]') == 3
				% Revolute around Z-axis
				Gm(rows,colsI) = GmI(3,:);
				Gm(rows,colsP) = GmP(3,:);
			end
			if this.vel
				gmdot(rows) = this.qdot; % Y-axis rotation
			else
				error('only works with velocity-level Euler');
			end
		end
	end
end
