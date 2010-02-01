#Manifold Class

setClass("manifold",
	representation(	name="character",	
					atlas="matrix",
					distance="function",
					origin="matrix"
	)
)
