// custom javascript
$(document).ready(() => {
  console.log("document ready!");
});
const csrftoken = getCookie('csrftoken');
// Function to call the function run_task in Views.py
$(".button").on("click", function () {
  $("#preLoader").fadeIn();
  console.log("button clicked");
  const filepath = $("#filepath").text() ? $("#filepath").text() : "none";
  console.log(filepath);
  $.ajax({
    url: "/tasks/",
    data: { type: $(this).data("type"), filepath: filepath.toString() },
    method: "POST",
  })
    .done((res) => {
      console.log(res);
      if (jQuery.isEmptyObject(res)) {
        console.log("reloading.")
        window.location.reload("/import/")
      }
      const html = `
      <tr>
      <td>${res.task_num}</td>
      <td>${res.task_status}</td>
      <td>${res.task_result}</td>
      </tr>`;
      $("#tasks").append(html);
      if (res.task_status === "Form Errors") {
        $("#save_Proceed").prop("disabled", true);
        $("#save_Proceed").addClass("disabled");
        $("#next_to_confirm").toggleClass("visible");
      }
      if (res.task_status === "Form is Valid") {
        $("#next_to_confirm").toggleClass("visible");
      }
      const taskStatus = res.status;
      $("#preLoader").fadeOut();
    })
    .fail((err) => {
      console.log(err);
    });
});

$("#save_Proceed").on("click", function () {
  $("#preLoader").fadeIn();
  $.ajax({
    url: "/tasks/proceed",
    data: {},
    method: "POST",
    headers: { 'X-CSRFToken': csrftoken },
  })
    .done((res) => {
      $("#preLoader").fadeOut();
      save_data(res);
    })
    .fail((err) => {
      console.log(err);
    });
});

$("#stop_Proceed").on("click", function () {
  $.ajax({
    url: "/tasks/cancel",
    data: {},
    method: "POST",
  })
    .done((res) => {
      window.alert("task canceled");
      location.reload();
    })
    .fail((err) => {
      console.log(err);
    });

});

save_data = (res) => {
  $.ajax({
    url: "/tasks/cancel",
    data: {},
    method: "POST",
  })
    .done((res) => {
      window.alert("Data Saved!");
      location.reload();
    })
    .fail((err) => {
      console.log(err);
    });

};