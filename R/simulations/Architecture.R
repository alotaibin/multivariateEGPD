# Neural network that transforms data into d-dimensional summary statistics
juliaEval('
     import NeuralEstimators: ResidualBlock
     function initializenetwork(n::Integer, d::Integer; residuals::Bool = true)
     
        # Expert summary statistics
        S(Z) = log.(samplesize(Z)) .- 7.7f0
     
        # Width of each hidden layer
        w = 128 
        
        # Input dimension to phi 
        w_phi = w + !isnothing(S)
        
        if residuals 
            # Architecture using residual connections and layer normalization
            function ResidualBlock(dim::Integer)
              # An MLP residual block: Dense -> Norm -> ReLU -> Add input
              layer = Chain(Dense(dim, dim), LayerNorm(dim), relu)
              return SkipConnection(layer, +)
            end
            psi = Chain(
                Dense(n, w), LayerNorm(w), relu,
                ResidualBlock(w),
                ResidualBlock(w)
            )
    
            phi = Chain(
                ResidualBlock(w_phi),
                ResidualBlock(w_phi),
                LayerNorm(w_phi), relu,
                Dense(w_phi, d)
            )
        else 
            # Simple architecture
            psi = Chain(Dense(n, w, relu), Dense(w, w, relu), Dense(w, w, relu))
            phi = Chain(Dense(w_phi, w, relu), Dense(w, w, relu), Dense(w, d))
        end 

        network = DeepSet(psi, phi; S = S)
        
        return network
    end 
')

n <- 2L # dimension of each data replicate (bivariate)
d <- 6L # dimension of the parameter vector 
J <- 5L # number of ensemble components for the NBE

# Initialize NPE 
NPE <- juliaLet('

  # Number of parameters in the model

  # Distribution used to approximate the posterior
  q = NormalisingFlow(d, 2d) 
  
  # Neural network that transforms data into summary statistics
  network = initializenetwork(n, 2d)

  # Neural posterior estimator
  NPE = PosteriorEstimator(q, network)
  ', d = d, n = n)

# Initialize NBE
NBE <- juliaLet('
    NBE = Ensemble(
    [PointEstimator(initializenetwork(n, d)) for j in 1:J]
    )  
    ', n = n, d = d, J = J)
