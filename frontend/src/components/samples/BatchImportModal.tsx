/**
 * Batch Import Modal - Import multiple samples from CSV
 */
import type React from 'react'
import { useState } from 'react'
import { Modal, Upload, Button, Steps, Table, Alert, Typography, Space, Tag, Progress, Card, message } from 'antd'
import {
  InboxOutlined,
  CloudUploadOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  DownloadOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons'
import type { UploadFile, UploadProps } from 'antd'
import sampleService from '@/services/sample.service'
import { toast } from '@/utils/notification'

const { Dragger } = Upload
const { Title, Text, Paragraph } = Typography

interface BatchImportModalProps {
  open: boolean
  projectId: string
  onClose: () => void
  onSuccess: () => void
}

interface ParsedSample {
  sample_id: string
  run_id?: string
  group_name?: string
  layout?: string
  batch_id?: string
  status?: 'pending' | 'success' | 'error'
  error?: string
}

export const BatchImportModal: React.FC<BatchImportModalProps> = ({ open, projectId, onClose, onSuccess }) => {
  const [currentStep, setCurrentStep] = useState(0)
  const [fileList, setFileList] = useState<UploadFile[]>([])
  const [parsedSamples, setParsedSamples] = useState<ParsedSample[]>([])
  const [uploading, setUploading] = useState(false)
  const [uploadProgress, setUploadProgress] = useState(0)

  const handleReset = () => {
    setCurrentStep(0)
    setFileList([])
    setParsedSamples([])
    setUploadProgress(0)
  }

  const handleClose = () => {
    handleReset()
    onClose()
  }

  // Parse CSV file
  const parseCSV = (file: File): Promise<ParsedSample[]> => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()

      reader.onload = (e) => {
        try {
          const text = e.target?.result as string
          const lines = text.split('\n').filter((line) => line.trim())

          if (lines.length < 2) {
            reject(new Error('CSV file must contain headers and at least one data row'))
            return
          }

          // Parse headers
          const headers = lines[0].split(',').map((h) => h.trim().toLowerCase())
          const sampleIdIndex = headers.indexOf('sample_id')

          if (sampleIdIndex === -1) {
            reject(new Error('CSV must contain a "sample_id" column'))
            return
          }

          // Parse data rows
          const samples: ParsedSample[] = lines.slice(1).map((line, index) => {
            const values = line.split(',').map((v) => v.trim())
            const sample: ParsedSample = {
              sample_id: values[sampleIdIndex] || `Sample_${index + 1}`,
              status: 'pending',
            }

            // Parse optional fields
            const runIdIndex = headers.indexOf('run_id')
            const groupIndex = headers.indexOf('group_name')
            const layoutIndex = headers.indexOf('layout')
            const batchIndex = headers.indexOf('batch_id')

            if (runIdIndex !== -1) {
              sample.run_id = values[runIdIndex] || undefined
            }
            if (groupIndex !== -1) {
              sample.group_name = values[groupIndex] || undefined
            }
            if (layoutIndex !== -1) {
              sample.layout = values[layoutIndex] || undefined
            }
            if (batchIndex !== -1) {
              sample.batch_id = values[batchIndex] || undefined
            }

            return sample
          })

          resolve(samples.filter((s) => s.sample_id))
        } catch (error: any) {
          reject(new Error(`Failed to parse CSV: ${error.message}`))
        }
      }

      reader.onerror = () => {
        reject(new Error('Failed to read file'))
      }

      reader.readAsText(file)
    })
  }

  const uploadProps: UploadProps = {
    name: 'file',
    multiple: false,
    fileList,
    maxCount: 1,
    accept: '.csv',
    beforeUpload: async (file) => {
      try {
        // Parse the CSV file
        const samples = await parseCSV(file)
        setParsedSamples(samples)
        setFileList([file as UploadFile])
        setCurrentStep(1)
        message.success(`Parsed ${samples.length} samples from CSV`)
      } catch (error: any) {
        message.error(error.message)
      }
      return false // Prevent auto upload
    },
    onRemove: () => {
      setFileList([])
      setParsedSamples([])
      setCurrentStep(0)
    },
  }

  const handleImport = async () => {
    if (parsedSamples.length === 0) {
      message.warning('No samples to import')
      return
    }

    setUploading(true)
    setCurrentStep(2)
    setUploadProgress(0)

    try {
      const formData = new FormData()
      formData.append('project_id', projectId)

      // Create CSV content from parsed samples
      const csvContent = [
        'sample_id,run_id,group_name,layout,batch_id',
        ...parsedSamples.map(
          (s) => `${s.sample_id},${s.run_id || ''},${s.group_name || ''},${s.layout || ''},${s.batch_id || ''}`,
        ),
      ].join('\n')

      const blob = new Blob([csvContent], { type: 'text/csv' })
      formData.append('file', blob, 'samples.csv')

      // Simulate progress
      const progressInterval = setInterval(() => {
        setUploadProgress((prev) => Math.min(prev + 10, 90))
      }, 200)

      const result = await sampleService.importFromCSV(projectId, fileList[0] as unknown as File)

      clearInterval(progressInterval)
      setUploadProgress(100)

      // Update sample statuses
      const updatedSamples = parsedSamples.map((s) => ({
        ...s,
        status: 'success' as const,
      }))
      setParsedSamples(updatedSamples)

      setTimeout(() => {
        toast.success(`Successfully imported ${result.count || parsedSamples.length} samples`)
        onSuccess()
        handleClose()
      }, 1000)
    } catch (error: any) {
      setUploadProgress(0)

      // Mark all as error
      const updatedSamples = parsedSamples.map((s) => ({
        ...s,
        status: 'error' as const,
        error: error.message,
      }))
      setParsedSamples(updatedSamples)

      toast.error(`Import failed: ${error.message}`)
      setCurrentStep(1) // Go back to preview
    } finally {
      setUploading(false)
    }
  }

  const downloadTemplate = () => {
    const template = `sample_id,run_id,group_name,layout,batch_id
Sample001,Run001,Control,PE,Batch1
Sample002,Run001,Treatment,PE,Batch1
Sample003,Run002,Control,SE,Batch2`

    const blob = new Blob([template], { type: 'text/csv' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = 'sample_import_template.csv'
    a.click()
    URL.revokeObjectURL(url)

    message.success('Template downloaded')
  }

  const columns = [
    {
      title: 'Sample ID',
      dataIndex: 'sample_id',
      key: 'sample_id',
      width: 150,
      render: (text: string) => <Text strong>{text}</Text>,
    },
    {
      title: 'Run ID',
      dataIndex: 'run_id',
      key: 'run_id',
      width: 120,
      render: (text?: string) => text || <Text type="secondary">-</Text>,
    },
    {
      title: 'Group',
      dataIndex: 'group_name',
      key: 'group_name',
      width: 120,
      render: (text?: string) => (text ? <Tag>{text}</Tag> : <Text type="secondary">-</Text>),
    },
    {
      title: 'Layout',
      dataIndex: 'layout',
      key: 'layout',
      width: 80,
      render: (text?: string) => (text ? <Tag color="blue">{text}</Tag> : <Text type="secondary">-</Text>),
    },
    {
      title: 'Batch ID',
      dataIndex: 'batch_id',
      key: 'batch_id',
      width: 100,
      render: (text?: string) => text || <Text type="secondary">-</Text>,
    },
    {
      title: 'Status',
      dataIndex: 'status',
      key: 'status',
      width: 100,
      render: (status?: string, _record?: ParsedSample) => {
        if (status === 'success') {
          return (
            <Tag icon={<CheckCircleOutlined />} color="success">
              Success
            </Tag>
          )
        }
        if (status === 'error') {
          return (
            <Tag icon={<CloseCircleOutlined />} color="error">
              Failed
            </Tag>
          )
        }
        return <Tag color="default">Pending</Tag>
      },
    },
  ]

  return (
    <Modal title="Batch Import Samples" open={open} onCancel={handleClose} width={900} footer={null} destroyOnClose>
      <Steps
        current={currentStep}
        items={[
          { title: 'Upload CSV', icon: <CloudUploadOutlined /> },
          { title: 'Preview', icon: <InfoCircleOutlined /> },
          { title: 'Import', icon: <CheckCircleOutlined /> },
        ]}
        style={{ marginBottom: 24 }}
      />

      {/* Step 0: Upload */}
      {currentStep === 0 && (
        <Space direction="vertical" size="large" style={{ width: '100%' }}>
          <Alert
            message="CSV Format Requirements"
            description={
              <div>
                <Paragraph style={{ marginBottom: 8 }}>
                  Your CSV file must include a <Text code>sample_id</Text> column. Optional columns:
                </Paragraph>
                <ul style={{ marginBottom: 8 }}>
                  <li>
                    <Text code>run_id</Text> - Sequencing run identifier
                  </li>
                  <li>
                    <Text code>group_name</Text> - Sample group (e.g., Control, Treatment)
                  </li>
                  <li>
                    <Text code>layout</Text> - Sequencing layout (PE or SE)
                  </li>
                  <li>
                    <Text code>batch_id</Text> - Batch identifier
                  </li>
                </ul>
                <Button type="link" icon={<DownloadOutlined />} onClick={downloadTemplate} style={{ padding: 0 }}>
                  Download Template CSV
                </Button>
              </div>
            }
            type="info"
            showIcon
          />

          <Dragger {...uploadProps}>
            <p className="ant-upload-drag-icon">
              <InboxOutlined style={{ fontSize: 48, color: '#1890ff' }} />
            </p>
            <p className="ant-upload-text">Click or drag CSV file to this area to upload</p>
            <p className="ant-upload-hint">
              Upload a CSV file containing sample information. The file will be parsed automatically.
            </p>
          </Dragger>
        </Space>
      )}

      {/* Step 1: Preview */}
      {currentStep === 1 && (
        <Space direction="vertical" size="large" style={{ width: '100%' }}>
          <Alert
            message={`Found ${parsedSamples.length} samples`}
            description="Please review the samples below before importing."
            type="success"
            showIcon
          />

          <Table
            columns={columns}
            dataSource={parsedSamples}
            rowKey="sample_id"
            pagination={{ pageSize: 10 }}
            scroll={{ y: 400 }}
            size="small"
          />

          <Space style={{ width: '100%', justifyContent: 'flex-end' }}>
            <Button onClick={handleReset}>Cancel</Button>
            <Button type="primary" onClick={handleImport} loading={uploading}>
              Import {parsedSamples.length} Samples
            </Button>
          </Space>
        </Space>
      )}

      {/* Step 2: Importing */}
      {currentStep === 2 && (
        <Space direction="vertical" size="large" style={{ width: '100%' }}>
          <Card>
            <Space direction="vertical" size="large" style={{ width: '100%', textAlign: 'center' }}>
              <CloudUploadOutlined style={{ fontSize: 64, color: '#1890ff' }} />
              <Title level={4}>Importing Samples...</Title>
              <Progress percent={uploadProgress} status={uploadProgress === 100 ? 'success' : 'active'} />
              <Text type="secondary">
                {uploadProgress < 100
                  ? `Uploading ${parsedSamples.length} samples to the server...`
                  : 'Import completed successfully!'}
              </Text>
            </Space>
          </Card>
        </Space>
      )}
    </Modal>
  )
}

export default BatchImportModal
