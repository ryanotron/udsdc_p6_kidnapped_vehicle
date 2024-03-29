/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * [x] TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * [x] TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
    num_particles = 8;  // [x] TODO: Set the number of particles
    
    std::default_random_engine gen;
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_th(theta, std[2]);
    double wt_init = 1.0;
    
    for (int i = 0; i < num_particles; i++) {
        Particle pi;
        pi.id = i;
        pi.x = dist_x(gen);
        pi.y = dist_y(gen);
        pi.theta = dist_th(gen);
        pi.weight = wt_init;
        
        particles.push_back(pi);
        weights.push_back(wt_init);
    }
    
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * [x] TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    
    // setup distributions
    std::default_random_engine gen;
    std::normal_distribution<double> dist_x(0, std_pos[0]);
    std::normal_distribution<double> dist_y(0, std_pos[1]);
    std::normal_distribution<double> dist_th(0, std_pos[2]);
    
    double x0, y0, th0;
    double v_per_w = velocity/yaw_rate;
    double dtheta = yaw_rate*delta_t;
    
    // update particle state one-by-one
    for (int i = 0; i < num_particles; i++) {
        x0 = particles[i].x;
        y0 = particles[i].y;
        th0 = particles[i].theta;

        // if yaw_rate is too small, v_per_w is NaN, then it's game over
        // so we need a fallback
        if (yaw_rate < -1.0e-9 || yaw_rate > 1.0e-9) {
            particles[i].x = x0 + v_per_w * (sin(th0 + dtheta) - sin(th0)) + dist_x(gen);
            particles[i].y = y0 + v_per_w * (cos(th0) - cos(th0 + dtheta)) + dist_y(gen);
        }
        else {
            particles[i].x = x0 + velocity*cos(th0)*delta_t + dist_x(gen);
            particles[i].y = y0 + velocity*sin(th0)*delta_t + dist_y(gen);
        }
        particles[i].theta = th0 + dtheta + dist_th(gen);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    
    int obsize = observations.size();
    // for each particle
    for (int i = 0; i < num_particles; i++) {
        Particle pi = particles[i];    
        vector<int> associations (obsize);
        vector<double> sense_x (obsize);
        vector<double> sense_y (obsize);
        double p_obs = 1.0;
        
        // for each observation
        for (unsigned int j = 0; j < observations.size(); j++) {
            // transform to map frame
            LandmarkObs obsj = observations[j];
            double x_m = obsj.x*cos(pi.theta) - obsj.y*sin(pi.theta) + pi.x;
            double y_m = obsj.x*sin(pi.theta) + obsj.y*cos(pi.theta) + pi.y;
            
            // find closest landmark in map
            int closest_id = 0;
            double closest_d = 100*sensor_range;
            double closest_x, closest_y;
            
            for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++) {
                Map::single_landmark_s markk = map_landmarks.landmark_list[k];
                double dee = dist(x_m, y_m, markk.x_f, markk.y_f);
                if (dee < closest_d) {
                    closest_d = dee;
                    closest_id = markk.id_i;
                    closest_x = markk.x_f;
                    closest_y = markk.y_f;
                }
            }
            
            if (closest_id < 0) {
                closest_x = x_m + sensor_range;
                closest_y = y_m + sensor_range;
            }
            
            associations[j] = closest_id;
            sense_x[j] = closest_x;
            sense_y[j] = closest_y;            
            
            // find measurement likelihood with gaussian mu = 0, 
            // x = predicted-observed, and sigma = sensor stdev
            // then immediately multiply together measurement likelihood for all 
            // observations
            
            // since a multivariate gaussian with zero off-diagonal cov matrix 
            // elements is just multiplication of univariate gaussians, we only 
            // have multiplication terms within and across observations
            // so, p_obs is a product "accumulator" for the particle under 
            // consideration
            
            // gaussian is defined in helper_functions.h
            p_obs *= gaussian(x_m - closest_x, 0, std_landmark[0])*gaussian(y_m - closest_y, 0, std_landmark[1]);
        }
        
        // set new weight
        pi.weight = p_obs;
        pi.associations = associations;
        pi.sense_x = sense_x;
        pi.sense_y = sense_y;
        particles[i] = pi;
        weights[i] = p_obs;
    }
}

void ParticleFilter::resample() {
  /**
   * [x] TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    
    std::default_random_engine gen;
    std::discrete_distribution<int> resample_dist(weights.begin(), weights.end());
    
    vector<Particle> new_particles;
    for (int i = 0; i < num_particles; i++) {
        int idx = resample_dist(gen);
        new_particles.push_back(particles[idx]);
    }
    
    particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
