function initialPivotTable(url) {
  var savecsv = null;

  $(document).ready(function () {
    console.log("loading data vitek");

    $(".submit_data").click(function () {
      collectAndSendData();
    });
  });

  const collectAndSendData = () => {
    var selected_data = [];
    $("input:checkbox[name=type]:checked").each(function () {
      selected_data.push($(this).val().toString());
    });
    var card_barcode = $("input[name=card_barcode]").val();
    var value_str = $("[data-name=data_process_value] option:selected")
      .val()
      .toString();
    var data_function_str = $("[data-name=data_function_type] option:selected")
      .val()
      .toString();

    var index_value = [];
    $("#sortable3 li").each(function () {
      index_value.push($(this).data("name"));
    });
    var index_values_str = index_value.toString();

    var column_value = [];
    $("#sortable2 li").each(function () {
      column_value.push($(this).data("name"));
    });
    var column_value_str = column_value.toString();

    var data = {
      selected_data: selected_data,
      values: value_str,
      columns: column_value_str,
      index: index_values_str,
      card_barcode: card_barcode,
      functions: data_function_str,
    };
    console.log(data);
    sendToServer(url, data);
  };

  const csrftoken = getCookie("csrftoken");
  const sendToServer = (data) => {
    console.log("send to server");

    $.ajax({
      url: url,
      type: "POST",
      headers: { "X-CSRFToken": csrftoken },
      data: data,
    })
      .done((response) => {
        console.log(response);
        data = response["table_html"];
        savecsv = response["table_csv"];
        $("#pivotable").html("");
        $("#pivotable").append(data);
      })
      .fail((XMLHttpRequest, textStatus, errorThrown) => {
        console.log(XMLHttpRequest, textStatus, errorThrown);
      });
  };

  $("#download_pt_as_csv").click(function () {
    console.log("clicked!");
    saveData(savecsv, "pivottable.csv");
  });
}
