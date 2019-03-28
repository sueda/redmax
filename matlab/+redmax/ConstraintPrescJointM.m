classdef ConstraintPrescJointM < redmax.Constraint
	% ConstraintPrescJointM Prescribes body motion in maximal coords
	% Only works with revolute joints
	
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
			% Assumes revolute around Y-axis
			Gm(rows,colsI) = GmI(2,:);
			Gm(rows,colsP) = GmP(2,:);
			if this.vel
				gmdot(rows) = this.qdot; % Y-axis rotation
			else
				error('only works with velocity-level Euler');
			end
		end
	end
end
