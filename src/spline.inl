// Given a time between 0 and 1, evaluates a cubic polynomial with
// the given endpoint and tangent values at the beginning (0) and
// end (1) of the interval.  Optionally, one can request a derivative
// of the spline (0=no derivative, 1=first derivative, 2=2nd derivative).
template <class T>
inline T Spline<T>::cubicSplineUnitInterval(
    const T& position0, const T& position1, const T& tangent0,
    const T& tangent1, double normalizedTime, int derivative) {
  
    if (!(0 <= derivative && derivative <= 2))
        return T(); //return something instead of breaking

    double t = normalizedTime;
    double t2 = t * t;
    double t3 = t * t * t; 

    double h00, h10, h01, h11;

    switch (derivative) {
    case 0:
        h00 = 2 * t3 - 3 * t2 + 1;
        h10 = t3 - 2 * t2 + t;
        h01 = -2 * t3 + 3 * t2;
        h11 = t3 - t2;
        break;
    case 1: 
        h00 = 6 * t2 - 6 * t;
        h10 = 3 * t2 - 4 * t + 1;
        h01 = -6 * t2 + 6 * t;
        h11 = 3 * t2 - 2 * t;
        break;
    case 2:
        h00 = 12 * t - 6;
        h10 = 6 * t + 4;
        h01 = -12 * t + 6;
        h11 = 6 * t - 2;
        break;
    }
    
    T value = h00 * position0 + h10 * tangent0;
            + h01 * position1 + h11 * tangent1;

    return value;
}

// Returns a state interpolated between the values directly before and after the
// given time.
template <class T>
inline T Spline<T>::evaluate(double time, int derivative) {
  // TODO (Animation) Task 1b

    if (knots.size() < 1) {
        return T();
    }

    else if (knots.size() == 1) {

        if (derivative >= 1) 
            return T();

        else
            return knots.begin()->second;
    }

    else if (knots.upper_bound(time) == knots.begin()) {

        if (derivative >= 1)
            return 0;

        else
            return knots.begin()->second;
        

    }

    else if (knots.upper_bound(time) == knots.end()) {
        
        if (derivative >= 1)
            return 0;

        else
            return knots.lower_bound(time)->second;
        
    }

    else {
        typename std::map<double, T>::iterator k0, k1, k2, k3, temp;
        k2 = knots.upper_bound(time);
        k1 = knots.lower_bound(time);

        double t0, t3;
        T p0, p3;

        if (k1 == knots.begin()) {
           
            t0 = 2 * (k1->first) - (k2->first);
            p0 = 2 * (k1->second) - (k2->second);
          
        }

        else {
            k0 = k1;
            k0--;

            t0 = k0->first;
            p0 = k0->second;
        }

        temp = k2; 
        temp++;

        if (temp == knots.end()) {
            t3 = 2 * (k2->first) - (k1->first);
            p3 = 2 * (k2->second) - (k1->second);
        }

        else {
            k3 = temp;
            t3 = k3->first;
            p3 = k3->second;
        }
        
        T p1, p2;

        p1 = k1->second;
        p2 = k2->second;

        double t1, t2;

        t1 = k1->first;
        t2 = k2->first;

        double t = t2 - t1;

        T m1 = (p2 - p0) / (t2 - t0) * t;
        T m2 = (p3 - p1) / (t3 - t1) * t;

        time = (time - t1) / t;

        return cubicSplineUnitInterval(p1, p2, m1, m2, time, derivative);

    }
}

// Removes the knot closest to the given time,
//    within the given tolerance..
// returns true iff a knot was removed.
template <class T>
inline bool Spline<T>::removeKnot(double time, double tolerance) {
  // Empty maps have no knots.
  if (knots.size() < 1) {
    return false;
  }

  // Look up the first element > or = to time.
  typename std::map<double, T>::iterator t2_iter = knots.lower_bound(time);
  typename std::map<double, T>::iterator t1_iter;
  t1_iter = t2_iter;
  t1_iter--;

  if (t2_iter == knots.end()) {
    t2_iter = t1_iter;
  }

  // Handle tolerance bounds,
  // because we are working with floating point numbers.
  double t1 = (*t1_iter).first;
  double t2 = (*t2_iter).first;

  double d1 = fabs(t1 - time);
  double d2 = fabs(t2 - time);

  if (d1 < tolerance && d1 < d2) {
    knots.erase(t1_iter);
    return true;
  }

  if (d2 < tolerance && d2 < d1) {
    knots.erase(t2_iter);
    return t2;
  }

  return false;
}

// Sets the value of the spline at a given time (i.e., knot),
// creating a new knot at this time if necessary.
template <class T>
inline void Spline<T>::setValue(double time, T value) {
  knots[time] = value;
}

template <class T>
inline T Spline<T>::operator()(double time) {
  return evaluate(time);
}
