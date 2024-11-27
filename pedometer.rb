class Filter
  COEFFICIENTS_LOW_0_HZ = {
    alpha: [1, -1.979133761292768, 0.979521463540373],
    beta:  [0.000086384997973502, 0.000172769995947004, 0.000086384997973502]
  }
  
  COEFFICIENTS_LOW_5_HZ = {
    alpha: [1, -1.80898117793047, 0.827224480562408],
    beta:  [0.095465967120306, -0.172688631608676, 0.095465967120306]
  }
  
  COEFFICIENTS_HIGH_1_HZ = {
    alpha: [1, -1.905384612118461, 0.910092542787947],
    beta:  [0.953986986993339, -1.907503180919730, 0.953986986993339]
  }

  def self.low_0_hz(data)
    filter(data, COEFFICIENTS_LOW_0_HZ)
  end

  def self.low_5_hz(data)
    filter(data, COEFFICIENTS_LOW_5_HZ)
  end

  def self.high_1_hz(data)
    filter(data, COEFFICIENTS_HIGH_1_HZ)
  end

  private

  def self.filter(input_signal, coefficients)
    filtered_output = [0,0]
    (2..input_signal.length-1).each do |i|
      filtered_output << coefficients[:alpha][0] *
                      (input_signal[i]     * coefficients[:beta][0] +
                       input_signal[i-1]   * coefficients[:beta][1] +
                       input_signal[i-2]   * coefficients[:beta][2] -
                       filtered_output[i-1] * coefficients[:alpha][1] -
                       filtered_output[i-2] * coefficients[:alpha][2])
    end
    filtered_output
  end
end

class Parser
  attr_reader :parsed_data

  def self.run(data)
    parser = Parser.new(data)
    parser.parse
    parser
  end

  def initialize(data)
    @data = data
  end

  def parse
    @parsed_data = @data.to_s.split(';').map { |sample| sample.split('|') }
                   .map { |axes| axes.map { |coords| coords.split(',').map(&:to_f) } }
    validate_data_format
    process_acceleration_components
  end

  private

  def validate_data_format
    unless @parsed_data.map { |x| x.map(&:length).uniq }.uniq == [[3]]
      raise 'Bad Input. Ensure data is properly formatted with x,y,z coordinates'
    end
  end

  def process_acceleration_components
    if @parsed_data.first.count == 1
      filtered_accl = @parsed_data.map(&:flatten).transpose.map do |total_accl|
        gravity_component = Filter.low_0_hz(total_accl)
        user_component = total_accl.zip(gravity_component).map { |a, b| a - b }
        [user_component, gravity_component]
      end

      @parsed_data = @parsed_data.length.times.map do |i|
        user = filtered_accl.map(&:first).map { |elem| elem[i] }
        grav = filtered_accl.map(&:last).map { |elem| elem[i] }
        [user, grav]
      end
    end
  end
end

class Processor
  attr_reader :dot_product_data, :filtered_data

  def self.run(data)
    processor = Processor.new(data)
    processor.dot_product
    processor.filter
    processor
  end

  def initialize(data)
    @data = data
  end

  def dot_product
    @dot_product_data = @data.map do |sample|
      user_accel = sample[0]
      grav_accel = sample[1]
      dot_product = user_accel[0] * grav_accel[0] + 
                   user_accel[1] * grav_accel[1] + 
                   user_accel[2] * grav_accel[2]
      dot_product.abs * 3.0
    end
  end

  def filter
    @filtered_data = @dot_product_data
    @filtered_data = Filter.low_5_hz(@filtered_data)
  end
end

class User
  GENDER = ['male', 'female']
  MULTIPLIERS = {'female' => 0.413, 'male' => 0.415}
  AVERAGES = {'female' => 70.0, 'male' => 78.0}

  attr_reader :gender, :height, :stride

  def initialize(gender = nil, height = nil, stride = nil)
    @gender = gender.to_s.downcase unless gender.to_s.empty?
    @height = Float(height) unless height.to_s.empty?
    @stride = Float(stride) unless stride.to_s.empty?

    raise 'Invalid gender' if @gender && !GENDER.include?(@gender)
    raise 'Invalid height' if @height && (@height <= 0)
    raise 'Invalid stride' if @stride && (@stride <= 0)

    @stride ||= calculate_stride
  end

  private

  def calculate_stride
    if gender && height
      MULTIPLIERS[@gender] * height
    elsif height
      height * (MULTIPLIERS.values.reduce(:+) / MULTIPLIERS.size)
    elsif gender
      AVERAGES[gender]
    else
      AVERAGES.values.reduce(:+) / AVERAGES.size
    end
  end
end

class Analyzer
  attr_reader :steps, :distance

  def self.run(filtered_data, user, trial)
    analyzer = Analyzer.new(filtered_data, user, trial)
    analyzer.analyze
    analyzer
  end

  def initialize(filtered_data, user, trial)
    @filtered_data = filtered_data
    @user = user
    @trial = trial
    @steps = 0
    @distance = 0
  end

  def analyze
    threshold = calculate_threshold
    @steps = count_steps(threshold)
    @distance = calculate_distance
  end

  private

  def calculate_threshold
    rates = []
    @filtered_data.each_cons(2) { |a, b| rates << (b - a) }
    
    mean_rate = rates.sum / rates.length
    std_dev = Math.sqrt(rates.map { |r| (r - mean_rate) ** 2 }.sum / rates.length)
    threshold = mean_rate * 0.3
    
    puts "\nThreshold Calculation Details:"
    puts "Mean rate: #{mean_rate.round(3)}"
    puts "Standard deviation: #{std_dev.round(3)}"
    puts "Threshold: #{threshold.round(3)}"
    
    threshold
  end

  def count_steps(threshold)
    step_count = 0
    min_peak_distance = 2
    last_peak_index = -min_peak_distance
    
    puts "Filtered data points: #{@filtered_data.length}"
    puts "All filtered values: #{@filtered_data.map { |v| v.round(3) }}"
    puts "Using threshold: #{threshold.round(3)}"
    
    rates = []
    @filtered_data.each_cons(2) { |a, b| rates << (b - a) }
    
    puts "\nRates of change: #{rates.map { |r| r.round(3) }}"
    
    (1..rates.length-2).each do |i|
      current_rate = rates[i]
      prev_rate = rates[i-1]
      next_rate = rates[i+1]
      
      is_peak = (prev_rate > threshold || current_rate > threshold) &&
                current_rate > next_rate &&
                (i - last_peak_index) > min_peak_distance
      
      if is_peak && !duplicate_peak?(i, last_peak_index, rates)
        step_count += 1
        last_peak_index = i
        puts "Step detected at index #{i} with rates: #{prev_rate.round(3)} -> #{current_rate.round(3)} -> #{next_rate.round(3)}"
      end
    end
    
    if rates.last > threshold && (rates.length - 1 - last_peak_index) > min_peak_distance
      step_count += 1
      puts "Final step detected with rate: #{rates.last.round(3)}"
    end
    
    puts "\nTotal peaks found: #{step_count}"
    step_count
  end

  def duplicate_peak?(current_index, last_peak_index, rates)
    return false if last_peak_index < 0
    current_rate = rates[current_index]
    last_peak_rate = rates[last_peak_index]
    (current_rate - last_peak_rate).abs < 1.0 && (current_index - last_peak_index) < 3
  end

  def calculate_distance
    @steps * @user.stride
  end
end

class Pipeline
  attr_reader :data, :user, :trial, :parser, :processor, :analyzer

  def self.run(data, user, trial)
    pipeline = Pipeline.new(data, user, trial)
    pipeline.feed
    pipeline
  end

  def initialize(data, user, trial)
    @data = data
    @user = user
    @trial = trial
  end

  def feed
    @parser = Parser.run(@data)
    @processor = Processor.run(@parser.parsed_data)
    @analyzer = Analyzer.run(@processor.filtered_data, @user, @trial)
  end
end

