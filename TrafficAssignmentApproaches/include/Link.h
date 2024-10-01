#ifndef LINK_H
#define LINK_H

#include <vector>
#include <iostream>
#include <cmath>


namespace TrafficAssignment {
  template <class T>
  struct Link {
    const int init, term, type;
    const long double power, capacity, length, free_flow_time, b, speed, toll;
    T flow;

    Link(int init, int term, T capacity, T length,
      T free_flow_time, T b, T power, T speed, T toll, int type) :
      init(init), term(term), type(type), capacity(capacity), length(length),
      free_flow_time(free_flow_time), b(b), power(power), speed(speed), toll(toll), flow(0) { };

    ~Link() { }

    // Calculates the delay function value based on the given flow and link characteristics
    T Delay(T temp_flow = -1) {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return free_flow_time;
      }
      return free_flow_time + b * free_flow_time * std::pow(temp_flow, power) / std::pow(capacity, power);
    }

    // Calculates the integrated delay function value based on the given flow and link characteristics
    T DelayInteg(T temp_flow = -1) {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return free_flow_time * temp_flow;
      }
      return free_flow_time * (temp_flow + b * capacity * std::pow(temp_flow / capacity, power + 1) / (power + 1));
    }

    // Calculates the differentiated delay function value based on the given flow and link characteristics
    T DelayDer(T temp_flow = -1) {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return 0;
      }
      return free_flow_time * b * power * std::pow(temp_flow / capacity, power - 1) / capacity;
    }

    // Calculates the double differentiated delay function value based on the given flow and link characteristics
    T DelaySecondDer(T temp_flow = -1) {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0) {
        return 0;
      }
      return free_flow_time * b * power * (power - 1) * std::pow(temp_flow, power - 2) / std::pow(capacity, power);
    }

    // Calculates the sum of provideded links delays
    static T GetLinksDelay(std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      T ans = 0;
      for (auto now : links_list) {
        ans += links[now].Delay();
      }
      return ans;
    }

    // Calculates the sum of provideded links differentiated delays
    static T GetLinksDelayDer(std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      T ans = 0;
      for (auto now : links_list) {
        ans += links[now].DelayDer();
      }
      return ans;
    }

    // Calculates the sum of provideded links doubly differentiated delays
    static T GetLinksDelaySecondDer(std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      T ans = 0;
      for (auto now : links_list) {
        ans += links[now].DelaySecondDer();
      }
      return ans;
    }
    
    // Checks if every link from the vector provided has a zero capacity
    static bool CheckNonZeroLinksCapacity(std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      bool ans = false;
      for (auto now : links_list) {
        if (links[now].capacity != 0) {
          ans = true;
        }
      }
      return ans;
    }
  };
}

#endif LINK_H