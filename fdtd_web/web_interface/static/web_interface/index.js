
window.onpopstate = function(event) {
  console.log(event.state.section);
  showSection(event.state.section);
}



document.addEventListener('DOMContentLoaded', function() {
    
    
    // Use buttons to toggle between views
    document.querySelector('#all-simulations-btn').addEventListener('click', () => load_simulations(response));
 
   

  
  });



  
  
  function add_simulation() {

    // Adjust the buttons in navbar
    document.querySelector('#all-simulations-btn').setAttribute('class','nav-link');
    document.querySelector('#add-simulation-btn').setAttribute('class','nav-link active');
    document.querySelector('#about-btn').setAttribute('class','nav-link');

    // Show the view where all simulations are shown and hide other views
    //document.querySelector('#simulation-view').style.display = 'none';
    document.querySelector('#all-simulations-view').style.display = 'none';
    document.querySelector('#add-simulation-view').style.display = 'block';
    document.querySelector('#about-view').style.display = 'none';

    
    const section = this.dataset.section;
   
    history.pushState({section: section},"", `${section}`);
    showSection(section);


    var modelInput;
      //Get the value of the select...
      document.querySelector('#ModelInput').addEventListener("input",function(){

        modelInput = document.querySelector('#ModelInput').options[document.querySelector('#ModelInput').selectedIndex].value;
        //console.log(modelInput);

        //Get the div tag...
        var device_model = document.querySelector('#input-model-group');
        //console.log(device_model);
        device_model.innerHTML = ""; //Empty the div before adding them..

        if(modelInput == "csv") //Put a File field in the div..
        {
          

          //Create an input element to be input to the div class
          var fileInput = document.createElement('input');
          fileInput.setAttribute('class','form-control');
          fileInput.setAttribute('type','file');
          fileInput.setAttribute('id','sim_input_filepath');

          //Create a label for the file
          var label_fileInput = document.createElement('label');
          label_fileInput.setAttribute('for',fileInput.getAttribute('id'));
          label_fileInput.setAttribute('class','form-label');
          label_fileInput.innerHTML = "Place the input file (csv) in this part";
          //Append to the div class
          device_model.appendChild(label_fileInput);
          device_model.appendChild(fileInput);

        
        }
        else if(modelInput == "manual") //Put the Input group here for manual parameters
        {
          //Get the number of models
          var input_param_row = document.createElement('div');
          input_param_row.setAttribute('class','row');
          input_param_row.id = "input_param";
          
          //Create the fmax input
          var fmax_div = document.createElement('div');
          fmax_div.setAttribute('class','col-4');
          fmax_div.id= "fmax_div";
          var fmax = document.createElement('input');
          fmax.type = "number";
          fmax.min = 0;
          fmax.setAttribute('class','form-control');
          fmax.required = true;
          fmax.id = "fmax";
          fmax.placeholder = "Frequency in Hz...";

          var label_fmax = document.createElement('label');
          label_fmax.innerHTML = "Frequency of Interest (in Hz)";
          label_fmax.setAttribute('for',fmax.id);
          label_fmax.setAttribute('class','form-label');
          
          fmax_div.appendChild(label_fmax);

          fmax_div.appendChild(fmax);
    
          
          
          input_param_row.appendChild(fmax_div);

          //Create the source type input
          var source_type = document.createElement('select');
          source_type.setAttribute('class','form-select form-select-lg mb-3 required');
          source_type.id = "source_type";
          source_type.required = true;
          source_type.innerHTML = `
                                  <option selected></option>
                                  <option value="gaussian"> Gaussian Pulse Source</option>
                                  <option value="sine"> Sinusoidal Source</option>
                                  <option value="square"> Rectangular Pulse Source</option>
                                  <option value="modulatedsine"> Modulated Sine Pulse Source</option>
                                  `;

          var source_type_div = document.createElement('div');
          source_type_div.setAttribute('class','col-4');
          source_type_div.innerHTML = "<label>Source Type</label>";
          source_type_div.appendChild(source_type);
          input_param_row.appendChild(source_type_div);

          //Create n_model input
          var n_model = document.createElement('input');
          n_model.type = "number";
          n_model.min = 1;
          n_model.setAttribute('class','form-control');
          n_model.required = true;
          n_model.id = "n_model";

          var nmodel_div = document.createElement('div');
          nmodel_div.innerHTML = "<label>Number of Device Layer/s</label>";
          nmodel_div.setAttribute('class','col-3');

          nmodel_div.appendChild(n_model);
          input_param_row.appendChild(nmodel_div);

          device_model.appendChild(input_param_row);

          var layers_div = document.createElement('div');
          layers_div.setAttribute('class','col-10');

          n_model.addEventListener('input',() => {

              //Clear the html info before repeating...
              layers_div.innerHTML="";
              for(let i = 1; i <= n_model.value; i++)
              {
                //Create a row for each layer
                var layer_heading = document.createElement('h6');
                layer_heading.innerHTML = "Layer " + i +":";
                var row = document.createElement('div');
                row.setAttribute('class','row');
                row.id = "row-" + i;

                

                //Create layer size
                var layer_size = document.createElement('input');
                layer_size.type = "number";
                layer_size.min = 0;
                layer_size.setAttribute('class','form-control ');
                layer_size.required = true;
                layer_size.step = "any";
                layer_size.id = "layer-" + i;


                var layer_size_div = document.createElement('div');
                layer_size_div.setAttribute('class','form-floating col-4');
                var layer_size_label = document.createElement('label');
                layer_size_label.innerHTML = "Size of layer (in meters)";
                layer_size_div.appendChild(layer_size);
                layer_size_div.appendChild(layer_size_label);
                row.appendChild(layer_size_div);

                //Create magnetic permeability
                var mu = document.createElement('input');
                mu.type = "number";
                mu.min = 0;
                mu.setAttribute('class','form-control ');
                mu.required = true;
                mu.step = "any";
                mu.id = "mu-" + i;


                var mu_div = document.createElement('div');
                mu_div.setAttribute('class','form-floating col-4');
                var mu_label = document.createElement('label');
                mu_label.innerHTML = "Relative magnetic permeability (μ_r)";
                mu_div.appendChild(mu);
                mu_div.appendChild(mu_label);
                row.appendChild(mu_div);

                //Create electric permittivity
                var epsilon = document.createElement('input');
                epsilon.type = "number";
                epsilon.min = 0;
                epsilon.setAttribute('class','form-control ');
                epsilon.required = true;
                epsilon.step = "any";
                epsilon.id = "epsilon-" + i;


                var epsilon_div = document.createElement('div');
                epsilon_div.setAttribute('class','form-floating col-4');
                var epsilon_label = document.createElement('label');
                epsilon_label.innerHTML = "Relative electric permittivity (ε_r)";
                epsilon_div.appendChild(epsilon);
                epsilon_div.appendChild(epsilon_label);
                row.appendChild(epsilon_div);

                layers_div.appendChild(layer_heading);
                layers_div.appendChild(row);
              }

          });

          device_model.appendChild(layers_div);
        }
        else{ //When the blank is selected
          device_model = document.querySelector('#input-model-group');
          //console.log(device_model);
          device_model.innerHTML = ""; //Empty the div before adding them..
          var alert_file = document.createElement('div');
          alert_file.setAttribute('class','alert alert-danger');
          alert_file.setAttribute('role','alert');
          alert_file.innerHTML = "NOTE: Please enter a mode of input for the simulation!";
          device_model.appendChild(alert_file);
        
        }
      });
      //Array used for getting the number of subdomains
      var subdomains = [1,2,4,8,16,32,64];
      //Get the div for the sim param
      var sim_param_group = document.querySelector('#sim-param-group');
      document.querySelector('#num_subdomains').value = 0;
      //Disable multithreading and number of subdomains by defauls
      document.querySelector('#multithreading-swtich').disabled = true;
      document.querySelector('#num_subdomains').disabled = true;
      document.querySelector('#currVal').innerHTML = subdomains[document.querySelector('#num_subdomains').value];
      document.querySelector('#algorithm').addEventListener("change",() => {
        var algorithmSelected  = document.querySelector('#algorithm')
                                .options[document.querySelector('#algorithm').selectedIndex].value;
        //console.log(algorithmSelected);
        
        console.log("algorithmSelected: " + algorithmSelected);
        if(algorithmSelected == "fdtd-schwarz")
        {
          //document.querySelector('#multithreading-swtich').checked = false;
          document.querySelector('#multithreading-swtich').disabled = false;
          document.querySelector('#num_subdomains').disabled = false;
          
        }
        else{
          document.querySelector('#multithreading-swtich').checked = false;
          document.querySelector('#multithreading-swtich').disabled = true;
          document.querySelector('#num_subdomains').value = 0;
          document.querySelector('#num_subdomains').disabled = true;
          document.querySelector('#currVal').innerHTML = subdomains[document.querySelector('#num_subdomains').value];
          //console.log("currVal: " + document.querySelector('#currVal').innerHTML + " | num_subdomains: " + document.querySelector('#num_subdomains').value);
        }
      });

      document.querySelector('#num_subdomains').addEventListener("change",() => {
        
        document.querySelector('#currVal').innerHTML = subdomains[document.querySelector('#num_subdomains').value];
        //console.log(document.querySelector('#currVal').innerHTML);

      });
    
    // Fetch all the forms we want to apply custom Bootstrap validation styles to
    const forms = document.querySelectorAll('.needs-validation');


    // Loop over them and prevent submission
    
   

    
    //fetch a POST request when the form is submitted
    document.querySelector('#sim_form').addEventListener('submit',(event)=>{

      // Prevent non-AJAX POST request to occur so that no wrong POST request is done
      event.preventDefault();

      //console.log(document.querySelector('#ModelInput').value);
      
      
      var form_valid = true;
      Array.from(forms).forEach(form => {
      
        if (!form.checkValidity()) {
          event.preventDefault()
          event.stopPropagation()
          form_valid = false;
        }
  
  
        form.classList.add('was-validated')
      })
      

      //console.log(document.querySelector('#ModelInput').value);
      //Fetch a POST request depending on the manual input..
      var dropdown = document.getElementById('ModelInput');
      //console.log(dropdown.value)

      var input_type = dropdown.value
      
     
      
      if (form_valid == true)
      {
        if(input_type == "csv")
        {
          var num_subdom = 0;
          if(document.querySelector('#algorithm').value == "fdtd")
            num_subdom = 1;
          else
            num_subdom = subdomains[document.querySelector('#currVal').innerHTML];
          fetch('new',{
            method: 'POST',
            body: JSON.stringify({
              username: document.querySelector('#username').value,
              user_email: document.querySelector('#user_email').value,
              sim_description: document.querySelector('#sim_description').value,
              input_type: document.querySelector('#ModelInput').value,
              input_filepath: document.querySelector('#sim_input_filepath').value,
              boundary_cond: document.querySelector('#boundary_cond').value,
              source_excitation: document.querySelector('#source_excitation').value,
              custom_name: document.querySelector('#custom_name').value,
              output_type: document.querySelector('#output_type').value,
              algorithm: document.querySelector('#algorithm').value,
              multithreading: document.querySelector('#multithreading-swtich').checked,
              num_subdomain: num_subdom,
              creation_datetime:  new Date().today() + " @ " + new Date().timeNow()
              
            })
          })
          .then(response => response.json())
          .then(result => {console.log(result)});
        }
        else if(input_type == "manual")
        {
          
          var num_layers = document.querySelector('#n_model').value;
          var layer_sizes = [];
          var mus = [];
          var epsilons = [];
          var layer_str = "#layer-";
          var mu_str = "#mu-";
          var epsilon_str = "#epsilon-";
          for(let i=1; i<= num_layers; i++){
            layer_sizes.push(document.querySelector(layer_str + i).value);
            mus.push(document.querySelector(mu_str + i).value);
            epsilons.push(document.querySelector(epsilon_str + i).value);
          }
          var num_subdom = 0;
          if(document.querySelector('#algorithm').value == "fdtd")
            num_subdom = 1;
          else
            num_subdom = subdomains[document.querySelector('#currVal').innerHTML];
          //Putting the CSRF token in the headers is a MUST when doing a POST request
          let headers = new Headers();
          headers.append('X-CSRFToken', csrftoken);
          fetch('new',{
            method: 'POST',
            //This is needed to include cookes (which has the CSRF token) for the POST request to be accepted
            credentials: 'same-origin',
            headers: headers,
            body: JSON.stringify({
              username: document.querySelector('#username').value,
              //csrfmiddlewaretoken: csrftoken,
              user_email: document.querySelector('#user_email').value,
              sim_description: document.querySelector('#sim_description').value,
              input_type: document.querySelector('#ModelInput').value,
              fmax: document.querySelector('#fmax').value,
              source_type: document.querySelector('#source_type').value,
              n_model: document.querySelector('#n_model').value,
              layer_size: layer_sizes,
              mu: mus,
              epsilon: epsilons,
              boundary_cond: document.querySelector('#boundary_cond').value,
              source_excitation: document.querySelector('#source_excitation').value,
              custom_name: document.querySelector('#custom_name').value,
              output_type: document.querySelector('#output_type').value,
              algorithm: document.querySelector('#algorithm').value,
              multithreading: document.querySelector('#multithreading-swtich').checked,
              num_subdomain: num_subdom,
              creation_datetime:  new Date().today() + " @ " + new Date().timeNow()
            })
          })
          .then(response => response.json())
          .then(result => {console.log(result)});
          
        }
        else{
          console.log("No input values detected");
          fetch('new',{
            method: 'POST',
            body: JSON.stringify({
              has_error: true,
              error_message: "Error: Incorrect value."
            })
          })
          .then(response => response.json())
          .then(result => {console.log("Entered POST request....")});
  
        }
      }
      else{

        //Get the time today
        var d = new Date().toLocaleTimeString();
        document.querySelector('#failed_submit_toast').textContent= d;

        //Show a toast saying that form validation failed...
        const toastLiveExample = document.getElementById('failed_submit_toast')
        const toast = new bootstrap.Toast(toastLiveExample)

        toast.show()
    
      }

   

      return false;

    });

  
    
    console.log("End of post request...");
  };

 // For todays date;
Date.prototype.today = function () { 
  return ((this.getDate() < 10)?"0":"") + this.getDate() +"/"+(((this.getMonth()+1) < 10)?"0":"") + (this.getMonth()+1) +"/"+ this.getFullYear();
}

// For the time now
Date.prototype.timeNow = function () {
   return ((this.getHours() < 10)?"0":"") + this.getHours() +":"+ ((this.getMinutes() < 10)?"0":"") + this.getMinutes() +":"+ ((this.getSeconds() < 10)?"0":"") + this.getSeconds();
}



  function showSection(section) {
    if(section == 'home'){
      // Show the view where all simulations are shown and hide other views
      //document.querySelector('#simulation-view').style.display = 'none';
      document.querySelector('#all-simulations-view').style.display = 'block';
      document.querySelector('#add-simulation-view').style.display = 'none';
      document.querySelector('#about-view').style.display = 'none';
    }

  }