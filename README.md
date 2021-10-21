# SBML_to_graph
A few *rudimentary* functions to generate metabolite graphs from SBML models. Uses ReFramed and one CarveMe function.

Includes a test jupyter-notebook (and scripts) with files to try it with.

### `create_spare_medium`

This function takes a given SBML model and simulates its growth on a given medium. Then returns a list of metabolites that are the "spent media", that is, the original media with the model's secretions and without the nutrients that the model has fully consumed.

Specifically, the "consumed" metabolites are specified by the user with the "to_remove" argument. But the products are automatically identified by parsing the ReFramed simulation results and saving the products of exchange type reactions with positive flux (if the flux is positive, metabolites are being excreted).

In metabolic models, there are normally three different types of pre-defined boundary reactions, according to the [COBRA documentation](https://cobrapy.readthedocs.io/en/latest/building_model.html): exchange, demand, and sink reactions. An exchange reaction is "a reversible reaction that adds to or removes an extracellular metabolite from the extracellular compartment", so that's why we are interested in those.

### `metabolite_translator`

`my_draw()` helper function. Uses "bigg_models_metabolites.txt" to turn BIGG IDs into pretty names for the graph.

### `reaction_parser`

`my_draw()` helper function. Applies weights (fluxes) and selects reactions by flux and by reversibility or irreversibility.

### `my_draw`

Draws graphs.
