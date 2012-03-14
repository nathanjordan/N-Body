/*
 * particle.h
 *
 *  Created on: Mar 13, 2012
 *      Author: njordan
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_


struct particle {

	float id, x , y, vx, vy, m;

	particle( float id ,float x ,float y ,float vx,float vy,float  m ) {
		this->id = id;
		this->x = x;
		this->y = y;
		this->vx = vx;
		this->vy = vy;
		this->m = m;
		}

	particle() {

		}

	particle& operator=( particle& p ) {
		this->id = id;
		this->x = x;
		this->y = y;
		this->vx = vx;
		this->vy = vy;
		this->m = m;
		return p;
		}

	};

#endif /* PARTICLE_H_ */
