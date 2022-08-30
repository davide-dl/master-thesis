include("../src/Dionysos.jl")

using OpenAIGym

# import Plots

using StaticArrays
using ..Dionysos
# using BenchmarkTools
DI = Dionysos
UT = DI.Utils
DO = DI.Domain
ST = DI.System
DM = DI.Map #domain_map
SR = DI.Relation #systemrelation
C = DI.Control

env = GymEnv(:CartPole, :v0)
for i ∈ 1:20
  T = 0
  R = run_episode(env, RandomPolicy()) do (s, a, r, s′)
    render(env)
    T += 1
  end
  @info("Episode $i finished after $T steps. Total reward: $R")
end
close(env)


# gravity = 9.8
# masscart = 1.0
# masspole = 0.1
# total_mass = masspole + masscart
# pole_length = 0.5  # actually half the pole's length
# polemass_length = masspole * spole_length
# force_mag = 10.0
# tau = 0.02  # seconds between state updates
# kinematics_integrator = "euler"
#
# # Angle at which to fail the episode
# theta_threshold_radians = 12 * 2 * math.pi / 360
# x_threshold = 2.4
#
# # Angle limit set to 2 * theta_threshold_radians so failing observation
# # is still within bounds.
# high = np.array(
#   [
#       self.x_threshold * 2,
#       np.finfo(np.float32).max,
#       self.theta_threshold_radians * 2,
#       np.finfo(np.float32).max,
#   ],
#   dtype=np.float32,
# )
#
# self.action_space = spaces.Discrete(2)
# self.observation_space = spaces.Box(-high, high, dtype=np.float32)
#
# self.screen = None
# self.clock = None
# self.isopen = True
# self.state = None
#
# self.steps_beyond_done = None
#
#   def step(self, action):
#       err_msg = f"{action!r} ({type(action)}) invalid"
#       assert self.action_space.contains(action), err_msg
#       assert self.state is not None, "Call reset before using step method."
#       x, x_dot, theta, theta_dot = self.state
#       force = self.force_mag if action == 1 else -self.force_mag
#       costheta = math.cos(theta)
#       sintheta = math.sin(theta)
#
#       # For the interested reader:
#       # https://coneural.org/florian/papers/05_cart_pole.pdf
#       temp = (
#           force + self.polemass_length * theta_dot**2 * sintheta
#       ) / self.total_mass
#       thetaacc = (self.gravity * sintheta - costheta * temp) / (
#           self.length * (4.0 / 3.0 - self.masspole * costheta**2 / self.total_mass)
#       )
#       xacc = temp - self.polemass_length * thetaacc * costheta / self.total_mass
#
#       if self.kinematics_integrator == "euler":
#           x = x + self.tau * x_dot
#           x_dot = x_dot + self.tau * xacc
#           theta = theta + self.tau * theta_dot
#           theta_dot = theta_dot + self.tau * thetaacc
#       else:  # semi-implicit euler
#           x_dot = x_dot + self.tau * xacc
#           x = x + self.tau * x_dot
#           theta_dot = theta_dot + self.tau * thetaacc
#           theta = theta + self.tau * theta_dot
#
#       self.state = (x, x_dot, theta, theta_dot)
#
#       done = bool(
#           x < -self.x_threshold
#           or x > self.x_threshold
#           or theta < -self.theta_threshold_radians
#           or theta > self.theta_threshold_radians
#       )
#
#       if not done:
#           reward = 1.0
#       elif self.steps_beyond_done is None:
#           # Pole just fell!
#           self.steps_beyond_done = 0
#           reward = 1.0
#       else:
#           if self.steps_beyond_done == 0:
#               logger.warn(
#                   "You are calling 'step()' even though this "
#                   "environment has already returned done = True. You "
#                   "should always call 'reset()' once you receive 'done = "
#                   "True' -- any further steps are undefined behavior."
#               )
#           self.steps_beyond_done += 1
#           reward = 0.0
#
#       return np.array(self.state, dtype=np.float32), reward, done, {}
