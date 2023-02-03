// custom javascript
$(document).ready(() => {
  console.log("document ready!");
});
var output_files = $("#output");

$("#id_file_field").change(function () {
  output_files.empty();
  var files = $(this).prop("files");
  var htmls = ``;
  jQuery.each(files, function (i, val) {
    htmls += `<p>${i + 1}-${val.name}</p>`;
  });
  output_files.append(htmls);
});
const csrftoken = getCookie("csrftoken");
if ($("#filepath").text()) {
  $("#progressbar span:first-child").toggleClass("bg-success");
}
// get list of validate passed files and be able to save to DB
var validatepassedfile = ["0"];
// Function to call the function run_task in Views.py

$(".button").on("click", function () {
  if ($("input:checked").length < 1) {
    alert("haven't select a file!")
    return
  } else {
    if ($(this).data("type") === "Validation") {
      $("#preLoader").fadeIn();
    }
    var filepathlist = [];
    $("input[name=uploadedfiles_select]:checked").each(function () {
      filepathlist.push($(this).val().toString());
    });

    const datamodel = $("#datamodel").text() ? $("#datamodel").text() : "none";

    $.ajax({
      url: "/import-VITEK/",
      data: {
        type: $(this).data("type"),
        filepathlist: filepathlist,
        datamodel: datamodel.toString(),
      },
      method: "POST",
      headers: { "X-CSRFToken": csrftoken },
    })
      .done((res) => {
        console.log(res)
        if (res.status === "Delete") {
          console.log($(this).find("p"))
          if (res.systemErr) {

            $(this).find("div").html("<div id='alert' class='badge bg-danger text-wrap'>file not existed!</div>")
          } else {

            $("#delete_cancel_task a:first-child").addClass("blink");
            $('input[name=uploadedfiles_select]:checked').addClass("not-visible");
            $('input[name=uploadedfiles_select]:checked').parent().addClass("not-visible");
            $(this).find("div").html("<div id='alert' class='badge bg-primary text-wrap'>file deleted!</div>")

          }
        } else {
          var validateResult = JSON.parse(res.validate_result.replace(/'/g, '"'));
          var validateReport = JSON.parse(
            res.file_report
              .replace(/'/g, '"')
              .replace(/"{/g, "{")
              .replace(/}"/g, "}")
          );
          console.log(validateReport);
          // JSON.parse;
          var f_list = Object.keys(validateResult);
          console.log(f_list)
          // ------------
          for (let i = 0; i < f_list.length; i++) {
            var error_num = 0;
            var warning_num = 0;
            var ew_description = { Error: [], Warning: [] };
            validateReport[f_list[i]].forEach((el) => {
              if (el["Error"].length > 0) {
                error_num++;
                ew_description["Error"].push(el["Error"].toString());
              }
              if (el["Warning"].length > 0) {
                warning_num++;
                ew_description["Warning"].push(el["Warning"].toString());
              }
            });


            const tr = `
          <tr>
          <td>${f_list[i]}</td>
          <td>${validateResult[f_list[i]]}</td>
          <td>${error_num.toString()}</td>
          <td>${warning_num.toString()}</td>
          <td><div id="upload_report">${JSON.stringify(ew_description)}</div></td>
          </tr>`;
            $("#tasks").append(tr);

            if (error_num === 0) {
              $('input[name=uploadedfiles_select]:checked').parent().addClass('text-success')
              $("#progressbar span:nth-child(2)").toggleClass("bg-success");
              validatepassedfile.push(f_list[i].toString().toLowerCase())
              // $("input[value=" + f_list[i] + "]").addClass("validatedfile_db")
              $(".confirmButton").prop("disabled", false);
            } else {
              $(".confirmButton").prop("disabled", true);
            }
          }
          // --------------
          $("#tasksreport").append(`***
        ${res.file_report} *** <br>`);

          $("#Import_step3").addClass("visible");

          $("#preLoader").fadeOut();
        }
      })
      .fail((XMLHttpRequest, textStatus, errorThrown) => {
        console.log(XMLHttpRequest, textStatus, errorThrown);
      });

  }
});
$("input[type=radio]").click(() => {
  console.log(validatepassedfile)
  console.log($('input[name=uploadedfiles_select]:checked').val().toString().toLowerCase())
  if (jQuery.inArray($('input[name=uploadedfiles_select]:checked').val().toString().toLowerCase(), validatepassedfile) > -1) {
    $(".confirmButton").prop("disabled", false);
  }
  else {
    $(".confirmButton").prop("disabled", true);
  }
})

$(".confirmButton").on("click", function () {
  if ($("input:checked").length < 1) {
    alert("haven't select a file!")
    return
  } {

    $("#preLoader").fadeIn();
    var filepathlist = [];
    $("input[name=uploadedfiles_select]:checked").each(function () {
      filepathlist.push($(this).val().toString());
    });
    const datamodel = $("#datamodel").text() ? $("#datamodel").text() : "none";
    $.ajax({
      url: "/import-VITEK/",
      method: "POST",
      data: {
        type: $(this).data("type"),
        filepathlist: filepathlist,
        datamodel: datamodel.toString(),
      },
      headers: { "X-CSRFToken": csrftoken },
    })
      .done((res) => {
        $("#preLoader").fadeOut();
        console.log(res);
        // if (res.task_status === "Form is Valid") {
        $("#Import_step4").addClass("visible");
        $("#progressbar span:nth-child(3)").toggleClass("bg-success");
        // }
        // const html = `<p>${res.status}</p>`;
        $("#mesg_save_Proceed").append(`<li>${res.status} saved!</li>`);
        // save_data(res);
      })
      .fail((err) => {
        console.log(err);
      });
  }
});
