// custom javascript
$(document).ready(() => {
  console.log("document ready!");
});
// ---------validating and uploading---------------------------------------------------//
var output_files = $("#output");
//Listing files before uploading
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

//---------Ajax call for validating, deleting, saving to DB--------------------------//
// Function to call Importhandler_VITEK in Views.py
//Button Event: could be validating-objects, deleting, validating-save-objects------//
$(".button").on("click", function () {
  // prevent from clicking without choosing a file
  if ($("input:checked").length < 1 || $("input:checked").hasClass('not-visible')) {
    alert("haven't select a file")
    return
  } else {
    if ($(this).data("type") === "Validation") {
      $("#preLoader").fadeIn();
    }
    // collecting input information: selected files...
    var select_file_list = [];
    $("input[name=uploadedfiles_select]:checked").each(function () {
      select_file_list.push($(this).val().toString());
    });

    const datamodel = $("#datamodel").text() ? $("#datamodel").text() : "none";
    // sending input data to server: files, process name, 
    $.ajax({
      url: "/import-VITEK/",
      data: {
        type: $(this).data("type"),
        select_file_list: select_file_list,
        datamodel: datamodel.toString(),
      },
      method: "POST",
      headers: { "X-CSRFToken": csrftoken },
    })
      .done((res) => {
        console.log(res)
        // process received DATA
        // ---------IN delete Case-----------
        if (res.status === "Delete") {
          if (res.systemErr) {
            $(this).find("div").html("<div id='alert' class='badge bg-danger text-wrap'>file not existed!</div>")
          } else {
            $("#delete_cancel_task a:first-child").addClass("blink");
            $('input[name=uploadedfiles_select]:checked').addClass("not-visible");
            $('input[name=uploadedfiles_select]:checked').parent().addClass("not-visible");
            $('input[name=uploadedfiles_select]:checked').prop("checked", false)
            $(this).find("div").html("<div id='alert' class='badge bg-primary text-wrap'>file deleted!</div>")
          }
        }

        else {
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
          // ------------Loop and parse each file's result mapping to the report table ----//
          var error_num = 0;
          var warning_num = 0;
          var ew_description = { Error: [], Warning: [] };
          for (let i = 0; i < f_list.length; i++) {
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
            // ----------setting file status: Case validating-save-objects without Errors ---------------
            if (res.status === "SavetoDB" && error_num === 0) {
              $.each($('input[name=uploadedfiles_select]:checked'), function (index, value) {

                if ($(this).val() === f_list[i]) {
                  if ($(this).parent().hasClass('text-danger')) {
                    $(this).parent().removeClass('text-danger')
                  }
                  $(this).parent().addClass('text-success')
                }
              })

              $("#preLoader").fadeOut();
              console.log(res)
              if (!$("#Import_step4").hasClass("visible")) {
                $("#Import_step4").addClass("visible");
              }
              $("#progressbar span:nth-child(3)").toggleClass("bg-success");
              $("#mesg_save_Proceed").append(`<li>${res.savefile} saved!</li>`);
            }
            //-----------setting file status: Case validating-objects without Error--------- 
            else if (res.status === "validating" && error_num === 0) {
              $.each($('input[name=uploadedfiles_select]:checked'), function (index, value) {

                if ($(this).val() === f_list[i]) {

                  if ($(this).parent().hasClass('text-danger')) {
                    $(this).parent().removeClass('text-danger')
                  }
                  $(this).parent().addClass('text-success')
                }
              })
              // if ($('input[value=' + f_list[i] + ']').parent().hasClass('text-danger')) {
              //   $('input[value=' + f_list[i] + ']').parent().removeClass('text-danger')
              // }
              // $('input[value=' + f_list[i] + ']').parent().addClass('text-success')
              $("#progressbar span:nth-child(2)").toggleClass("bg-success");
              validatepassedfile.push(f_list[i].toString().toLowerCase())
              // $("#confirmButton").prop("disabled", false);
              const tr = `
              <tr>
              <td>${f_list[i]}</td>
              <td>${validateResult[f_list[i]]}</td>
              <td>${error_num.toString()}</td>
              <td>${warning_num.toString()}</td>
              <td><div id="upload_report">${JSON.stringify(ew_description)}</div></td>
              </tr>`;
              $("#tasks").append(tr);
              $("#tasksreport").append(`***
              ${res.file_report} *** <br>`);


            }
            //-----------setting file status:  Case validating-objects or validating-save-objects with Error occurs 
            else {
              // $("#confirmButton").prop("disabled", true);
              $.each($('input[name=uploadedfiles_select]:checked'), function (index, value) {

                if ($(this).val() === f_list[i]) {

                  if ($(this).parent().hasClass('text-success')) {
                    $(this).parent().removeClass('text-success')
                  }
                  $(this).parent().addClass('text-danger')
                }
              })
              const tr = `
              <tr>
              <td>${f_list[i]}</td>
              <td>${validateResult[f_list[i]]}</td>
              <td>${error_num.toString()}</td>
              <td>${warning_num.toString()}</td>
              <td><div id="upload_report">${JSON.stringify(ew_description)}</div></td>
              </tr>`;
              $("#tasks").append(tr);
              $("#tasksreport").append(`***
              ${res.file_report} *** <br>`);
            }

            //--------Single File process fisnish-------------------------------------------//
          }
          //----------Files Looping End----------------------------------------------------//

          // --------------Case After click Validating
          if (res.status === "validating") {
            if (!$("#Import_step3").hasClass("visible")) {

              $("#Import_step3").addClass("visible");
            }
            if (error_num === 0) {
              $("#confirmButton").prop("disabled", false);
            }
          }
          // ------end preloader displaying 
          $("#preLoader").fadeOut();
        }
      })
      // If Errors when Ajax Call
      .fail((XMLHttpRequest, textStatus, errorThrown) => {
        console.log(XMLHttpRequest, textStatus, errorThrown);
      });

  }
});
//------Button Click Event End-------------------------------------------------------------------------//
// -----save button disable when choose files contains error. Meaning Only passing Validtion files can be save. 
$("input[type=checkbox]").click(() => {
  $("#confirmButton").prop("disabled", false);
  $.each($('input[name=uploadedfiles_select]:checked'), function (index, value) {
    if ($(this).parent().hasClass('text-danger')) {

      $("#confirmButton").prop("disabled", true);
    }
  })

})
