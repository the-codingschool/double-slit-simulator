class waveFunction():

  """
  Calculates the 2D wavefunction of an electron in the double slit experiment.


  © 2024 The Coding School, All rights reserved.
  """

  d: float
  """The distance between the two slits. The larger this distance gets, the more the interference pattern starts to dissapear as the peaks blend together and we recover more classical results."""

  distance_to_screen: float
  """The distance between the slits and the screen where the particles will be measured. The larger this distance gets, the more spread out the interference pattern gets."""

  measure_slit: bool
  """Whether to measure which slit the particle goes through or not. Measuring at the slits collapses the wavefunction there instead of at the screen, thereby destroying the wave like nature of the particles."""


  def __init__(self,d, distance_to_screen,measure_slit):

    """
    Initializes the double slit wavefunction.

    ###Parameters:
    <ul>
      <li>`d` (float): The distance between the two slits.</li>
      <li>`distance_to_screen` (float): The distance between the slits and the screen where the particles will be measured.</li>
      <li>`measure_slit` (bool): Whether to measure which slit the particle goes through or not.</li>
    </ul>
    """

    self.d = d
    self.distance_to_screen = distance_to_screen
    self.measure_slit= measure_slit

    if not self.measure_slit:
      self.values = np.linspace(-10,10,num=1000)
      self.norm = scipy.integrate.trapezoid(self.evaluate_unnormalized(self.values),self.values)
      self.probs = [0]
      self.probs.extend([scipy.integrate.trapezoid(self.evaluate(np.linspace(self.values[i],self.values[i+1],num=100)),np.linspace(self.values[i],self.values[i+1],num=100)) for i in range(999)])
      self.probs = self.probs/sum(self.probs)
    else:
      self.values = [-1*self.d/2,self.d/2]
      self.norm = 1
      self.probs = [0.5,0.5]

  def evaluate(self,x):
    """Returns the wavefunction probability distribution evaluated at a specific point. Only used for non-collapsed wavefunctions.

    ###Parameters:
    <ul>
      <li>`x` (float): The point to evaluate the wavefunction at.</li>
    </ul>
    
    ###Returns:
    The wavefunction evaluated at the provided point.


    ### NOTE:
    This evaluation method is not used to extract probabilities due to the nature of the
    probability distribution function.
    """


    if not self.measure_slit:
      return np.cos(np.pi * self.d* x/self.distance_to_screen)**2/self.norm
    else:
      if float(x)== float(-1*self.d/2):
        return 0.5
      elif float(x)== float(self.d/2):
        return 0.5
      else:
        return 0


  def evaluate_unnormalized(self,x):

    """Returns the unnormalized wavefunction probability distribution evaluated at a specific point.

    ###Parameters:
    <ul>
      <li>`x` (float): The point to evaluate the wavefunction at.</li>
    </ul>
    
    ###Returns:
    The unnormalized wavefunction evaluated at the provided point.


    ### NOTE:
    This evaluation method is only used to normalize the wavefunction.
    """

    return np.cos(np.pi * self.d * x/self.distance_to_screen)**2


  def measure(self):

    """Used to collapse the superposition and find the detected location of the electron.

    ###Parameters:
    None
    
    ###Returns:
    The x coordinate of the measured electron.
    """

    temp_value = np.random.choice(self.values, p=self.probs)
    if self.measure_slit:
      temp_value += np.random.normal(scale = 0.2)
    else:
      temp_value += np.random.uniform(low=-0.01,high=0.01)
    return temp_value



class doubleSlit():

  """
  Simulates the double slit experiment.

  The double slit experiment is a famous demonstration of quantum mechanics, showcasing the wave-particle duality of matter. \
  In this experiment, particles such as electrons are fired at a barrier with two narrow slits, and then detected on a screen positioned some distance behind the barrier. 
  The resulting pattern on the screen exhibits interference, which are characteristic of waves, despite the particles being fired one at a time.
  If one tries to measure which slit the particles travel through, the results on the screen drastically change to display clear particle like behavior,
  with no interference patterns. These results show the wave-particle duality of quantum physics.


  © 2024 The Coding School, All rights reserved.
  """

  slit_dist: float
  """The distance between the two slits. The larger this distance gets, the more the interference pattern starts to dissapear as the peaks blend together and we recover more classical results."""

  distance_to_screen: float
  """The distance between the slits and the screen where the particles will be measured. The larger this distance gets, the more spread out the interference pattern gets."""

  measure_slit: bool
  """Whether to measure which slit the particle goes through or not. Measuring at the slits collapses the wavefunction there instead of at the screen, thereby destroying the wave like nature of the particles."""

  screen_width: int
  """Number of bins along the x axis of the screen. Increases the resolution of
  the screen, but makes patterns harder to see unless the number of
  electrons is increased."""

  
  screen_height: int
  """Number of bins along the y axis of the screen. Increases the resolution of
  the screen, but makes patterns harder to see unless the number of
  electrons is increased."""

  detections_x: list
  """Where on the screen's x-axis each particle was measured. For internal use only. <b>SHOULD NOT BE MODIFIED BY USER</b>."""
  
  detections_y: list
  """Where on the screen's y-axis each particle was measured. For internal use only. <b>SHOULD NOT BE MODIFIED BY USER</b>."""
  
  wavefunction: waveFunction
  """The particle's wavefunction for a given experimental setup. For internal use only. <b>SHOULD NOT BE MODIFIED BY USER</b>."""



  def __init__(self,slit_dist = 1, distance_to_screen = 10, screen_width = 200, screen_height=100, measure_slit = False):

    """
    Initializes the double slit experiment.

    ###Parameters:
    <ul>
      <li>`slit_dist` (float, optional): The distance between the two slits. Defaults to `1`.</li>
      <li>`distance_to_screen` (float, optional): The distance between the slits and the screen where the particles will be measured. Defaults to `10`.</li>
      <li>`screen_width` (float, optional): Number of bins along the x axis of the screen. Increases the resolution of the screen, but makes patterns harder to see unless the number of electrons is increased. Defaults to `200`.</li>
      <li>`screen_height` (float, optional): Number of bins along the y axis of the screen. Increases the resolution of the screen, but makes patterns harder to see unless the number of electrons is increased. Defaults to `100`.</li>
      <li>`measure_slit` (bool, optional): Whether to measure which slit the particle goes through or not. Defaults to `False`.</li>
    </ul>
    """

    self.slit_dist = slit_dist
    self.distance_to_screen = distance_to_screen
    self.detections_x = []
    self.detections_y = []
    self.screen_width = screen_width
    self.screen_height = screen_height
    self.measure_slit = measure_slit

    self.wavefunction = waveFunction(self.slit_dist, self.distance_to_screen,self.measure_slit)


  def fire_electron(self):
    """Fires a single electron through the slits.  

    ###Parameters:
    None
    
    ###Returns:
    None

    ###Raises:
    `ValueError`
    Raised if any attributes of the `doubleSlit` object have been modified
    and no longer match those of the `waveFunction` object. See `clear_screen()`.
    """

    if self.slit_dist != self.wavefunction.d:
      raise ValueError("slit_dist attribute has been modified. Screen must be cleared.")
    elif self.distance_to_screen != self.wavefunction.distance_to_screen:
      raise ValueError("distance_to_screen attribute has been modified. Screen must be cleared.")

    detected_x = self.distance_to_screen*np.tan(self.wavefunction.measure())
    self.detections_x.append(self.wavefunction.measure())
    self.detections_y.append(np.random.normal(scale=1.7))


  def electron_beam(self, num_electrons = 5000):
    """Fires many electrons through the slits.  
    
    ###Parameters:
    <ul>
    <li>`num_electrons` (int): The number of electrons to fire. Defaults to 5000.</li>
    </ul>
    
    ###Returns:
    None

    ###Raises:
    `ValueError`
    Raised if any attributes of the `doubleSlit` object have been modified
    and no longer match those of the `waveFunction` object. See `clear_screen()`.
    """

    if self.slit_dist != self.wavefunction.d:
      raise ValueError("slit_dist attribute has been modified. Screen must be cleared.")
    elif self.distance_to_screen != self.wavefunction.distance_to_screen:
      raise ValueError("distance_to_screen attribute has been modified. Screen must be cleared.")
    for i in range(num_electrons):
      self.fire_electron()

 
  def show_screen(self):
    """Displays the screen. The range of the screen is [-10,10] along the x axis, and [-5,5] along the y axis.

    ###Parameters:
    None
    
    ###Returns:
    None
    """

    plt.hist2d(self.detections_x,self.detections_y,[self.screen_width,self.screen_height],range=[[-10,10],[-5,5]])
    plt.minorticks_on()
    plt.show()


  def clear_screen(self):
    """Clears the screen and recalculates the wavefunction. Useful if any attributes of the
    `doubleSlit` object are modified directly instead of creating a new object.
    
    ###Parameters:
    None
    
    ###Returns:
    None

    """

    self.detections_x = []
    self.detections_y = []
    self.wavefunction = waveFunction(self.slit_dist, self.distance_to_screen,self.measure_slit)


  def show_hist(self):
    """Plots a histogram of the x coordinates of the detected electrons.

    ###Parameters:
    None
    
    ###Returns:
    None

    """

    plt.hist(self.detections_x,bins=self.screen_width)
    plt.xlabel("Distance from center")
    plt.ylabel("Number of Electrons Detected")
    plt.minorticks_on()
    plt.show()
