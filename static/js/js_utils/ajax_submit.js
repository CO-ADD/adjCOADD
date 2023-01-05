$("#{some FormName}").on("submit", function (event) {
  event.preventDefault();
  create_entry();
});

// AJAX for posting
function create_entry() {
  console.log("sanity check!"); //
  $.ajax({
    url: "whereTopost", // the endpoint
    type: "POST", // http method
    data: { new_data: $("#{some FormName}").val() }, // data sent with the post request

    // handle a successful response
    success: function (json) {
      $("#{some FormName}").val(""); // remove the value from the input
      console.log(json); // log the returned json to the console
      console.log("success"); // another sanity check
    },

    // handle a non-successful response
    error: function (xhr, errmsg, err) {
      $("#form_result").html(
        "<div class='alert-box alert radius' data-alert>Oops! We have encountered an error: " +
          errmsg +
          " <a href='#' class='close'>&times;</a></div>"
      ); // add the error to the dom
      console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
    },
  });
}
