#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"
#include <cppad/cppad.hpp>
#include <chrono>
//#define DEBUG

// for convenience
using json = nlohmann::json;

bool is_initialized = false;

double cte_pred = 0;
double epsi_pred = 0;
double epsi = 0;
double delta = 0;
std::chrono::steady_clock::time_point last_time = std::chrono::steady_clock::now();

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);
  
  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }
  
  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }
  
  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

// convert miles per hour to meter per second
double convert_velocity(double speed_milesph) {
  double speed_meterps = 1609.34*speed_milesph/3600;
  return speed_meterps;
}

void transform_reference_system(const vector<double> ptsx_gloRF, const vector<double> ptsy_gloRF, double px, double py, double psi, vector<double> &refx_carRF, vector<double> &refy_carRF) {
  
  for (int i=0; i<ptsx_gloRF.size(); ++i) {
    refx_carRF.push_back((ptsx_gloRF[i]-px)*cos(psi)+(ptsy_gloRF[i]-py)*sin(psi));
    refy_carRF.push_back(-(ptsx_gloRF[i]-px)*sin(psi)+(ptsy_gloRF[i]-py)*cos(psi));
  }
}

int main() {
  uWS::Hub h;
  
  // MPC is initialized here!
  MPC mpc;
  
  
  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          
          std::chrono::steady_clock::time_point current_time = std::chrono::steady_clock::now();
          double elapsed_secs = std::chrono::duration_cast<std::chrono::microseconds>(current_time - last_time).count() /1000000.0;
          last_time = std::chrono::steady_clock::now();
          
          std::cout<<elapsed_secs<<std::endl;
          // j[1] is the data JSON object
          
          vector<double> refx_gloRF = j[1]["ptsx"];
          vector<double> refy_gloRF = j[1]["ptsy"];
          
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double v_ms = convert_velocity(v);
          
          vector<double> refx_carRF;
          vector<double> refy_carRF;
          
          transform_reference_system(refx_gloRF, refy_gloRF, px, py, psi, refx_carRF, refy_carRF);
          
          // converting to vector to eigen-vector
          double* ptrx = &refx_carRF[0];
          Eigen::Map<Eigen::VectorXd> posx_carRF_eigen(ptrx, refx_carRF.size());
          
          double* ptry = &refy_carRF[0];
          Eigen::Map<Eigen::VectorXd> posy_carRF_eigen(ptry, refy_carRF.size());
          
          auto coeffs = polyfit(posx_carRF_eigen, posy_carRF_eigen, 3);
          
          double x_eval = v_ms*elapsed_secs; // calculate trajectory for 100s to the future;
          double f = polyeval(coeffs, x_eval);
          // desired psi calculated as arctan of derivate of f
          double psi_des = CppAD::atan(3*x_eval*x_eval*coeffs[3]+2*x_eval*coeffs[2]+coeffs[1]);
          
          // in vehicle coordinates
          cte_pred = (f - 0) + (v_ms * CppAD::sin(epsi) * mpc.dt_);
          epsi_pred = (0 - psi_des) + v_ms * delta / mpc.Lf_ * mpc.dt_;
          
          Eigen::VectorXd state(6);
          state << x_eval, 0, 0, v, cte_pred, epsi_pred;
          
          auto vars = mpc.Solve(state, coeffs);
          epsi = vars[5];
          double steer_value = vars[6];
          double throttle_value = vars[7];
          
          json msgJson;
          msgJson["steering_angle"] = -steer_value;
          msgJson["throttle"] = throttle_value;
          
          //Display the MPC predicted trajectory
          msgJson["mpc_x"] = mpc.Xpred_;
          msgJson["mpc_y"] = mpc.Ypred_;
          
          //Display the reference line
          msgJson["next_x"] = refx_carRF;
          msgJson["next_y"] = refy_carRF;
          
          
          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          is_initialized = true;
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });
  
  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });
  
  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });
  
  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });
  
  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
