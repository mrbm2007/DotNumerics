namespace test2
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.button_fft = new System.Windows.Forms.Button();
            this.button_arma = new System.Windows.Forms.Button();
            this.button_ab = new System.Windows.Forms.Button();
            this.textBox_n = new System.Windows.Forms.TextBox();
            this.textBox_res = new System.Windows.Forms.TextBox();
            this.checkBox_arma = new System.Windows.Forms.CheckBox();
            this.button_clear = new System.Windows.Forms.Button();
            this.button1 = new System.Windows.Forms.Button();
            this.button2 = new System.Windows.Forms.Button();
            this.SuspendLayout();
            // 
            // button_fft
            // 
            this.button_fft.Location = new System.Drawing.Point(27, 41);
            this.button_fft.Name = "button_fft";
            this.button_fft.Size = new System.Drawing.Size(75, 23);
            this.button_fft.TabIndex = 0;
            this.button_fft.Text = "FFT";
            this.button_fft.UseVisualStyleBackColor = true;
            this.button_fft.Click += new System.EventHandler(this.button1_Click);
            // 
            // button_arma
            // 
            this.button_arma.Location = new System.Drawing.Point(27, 12);
            this.button_arma.Name = "button_arma";
            this.button_arma.Size = new System.Drawing.Size(75, 23);
            this.button_arma.TabIndex = 1;
            this.button_arma.Text = "ARMADILO";
            this.button_arma.UseVisualStyleBackColor = true;
            this.button_arma.Click += new System.EventHandler(this.button_arma_Click);
            // 
            // button_ab
            // 
            this.button_ab.Location = new System.Drawing.Point(27, 70);
            this.button_ab.Name = "button_ab";
            this.button_ab.Size = new System.Drawing.Size(75, 23);
            this.button_ab.TabIndex = 2;
            this.button_ab.Text = "A x B";
            this.button_ab.UseVisualStyleBackColor = true;
            this.button_ab.Click += new System.EventHandler(this.button3_Click);
            // 
            // textBox_n
            // 
            this.textBox_n.Location = new System.Drawing.Point(108, 73);
            this.textBox_n.Name = "textBox_n";
            this.textBox_n.Size = new System.Drawing.Size(44, 20);
            this.textBox_n.TabIndex = 3;
            this.textBox_n.Text = "3";
            // 
            // textBox_res
            // 
            this.textBox_res.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.textBox_res.Location = new System.Drawing.Point(12, 112);
            this.textBox_res.Multiline = true;
            this.textBox_res.Name = "textBox_res";
            this.textBox_res.ScrollBars = System.Windows.Forms.ScrollBars.Vertical;
            this.textBox_res.Size = new System.Drawing.Size(518, 309);
            this.textBox_res.TabIndex = 4;
            // 
            // checkBox_arma
            // 
            this.checkBox_arma.AutoSize = true;
            this.checkBox_arma.Checked = true;
            this.checkBox_arma.CheckState = System.Windows.Forms.CheckState.Checked;
            this.checkBox_arma.Location = new System.Drawing.Point(108, 16);
            this.checkBox_arma.Name = "checkBox_arma";
            this.checkBox_arma.Size = new System.Drawing.Size(66, 17);
            this.checkBox_arma.TabIndex = 5;
            this.checkBox_arma.Text = "Armadilo";
            this.checkBox_arma.UseVisualStyleBackColor = true;
            this.checkBox_arma.CheckedChanged += new System.EventHandler(this.checkBox1_CheckedChanged);
            // 
            // button_clear
            // 
            this.button_clear.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.button_clear.Location = new System.Drawing.Point(491, 399);
            this.button_clear.Name = "button_clear";
            this.button_clear.Size = new System.Drawing.Size(39, 22);
            this.button_clear.TabIndex = 6;
            this.button_clear.Text = "clear";
            this.button_clear.UseVisualStyleBackColor = true;
            this.button_clear.Click += new System.EventHandler(this.button_clear_Click);
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(174, 70);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 7;
            this.button1.Text = "Roots";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click_1);
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(255, 71);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(75, 23);
            this.button2.TabIndex = 8;
            this.button2.Text = "Eigen(A[])";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(542, 433);
            this.Controls.Add(this.button2);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.button_clear);
            this.Controls.Add(this.checkBox_arma);
            this.Controls.Add(this.textBox_res);
            this.Controls.Add(this.textBox_n);
            this.Controls.Add(this.button_ab);
            this.Controls.Add(this.button_arma);
            this.Controls.Add(this.button_fft);
            this.Name = "Form1";
            this.StartPosition = System.Windows.Forms.FormStartPosition.CenterScreen;
            this.Text = "DotNumeric Test";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button button_fft;
        private System.Windows.Forms.Button button_arma;
        private System.Windows.Forms.Button button_ab;
        private System.Windows.Forms.TextBox textBox_n;
        private System.Windows.Forms.TextBox textBox_res;
        private System.Windows.Forms.CheckBox checkBox_arma;
        private System.Windows.Forms.Button button_clear;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Button button2;
    }
}

